#include <gmsh.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <cmath>
#include <vector>
#include <string>

int main() {
  gmsh::initialize();
  gmsh::model::add("anim");

  gmsh::merge("Gengar.stl");
  gmsh::model::mesh::classifySurfaces(80 * M_PI / 180, true, true, M_PI);
  gmsh::model::mesh::createGeometry();

  std::vector<std::pair<int, int>> surfaces;
  gmsh::model::getEntities(surfaces, 2);

  std::vector<int> surfaceTags;
  for (auto s : surfaces) surfaceTags.push_back(s.second);

  gmsh::model::geo::synchronize();
  int loop = gmsh::model::geo::addSurfaceLoop(surfaceTags);
  gmsh::model::geo::addVolume({loop});
  gmsh::model::geo::synchronize();

  gmsh::option::setNumber("Mesh.MeshSizeMin", 1.4);
  gmsh::option::setNumber("Mesh.MeshSizeMax", 2.4);
  gmsh::model::mesh::generate(3);

  std::vector<std::size_t> nodeTags;
  std::vector<double> nodeCoords, tmp;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, tmp);

  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 3);

  int tetraIndex = -1;
  for (int i = 0; i < elementTypes.size(); i++)
    if (elementTypes[i] == 4) { tetraIndex = i; break; }

  int pointCount = nodeTags.size();

  std::size_t maxTag = 0;
  for (int i = 0; i < pointCount; i++)
    if (nodeTags[i] > maxTag) maxTag = nodeTags[i];

  std::vector<int> tagToIndex(maxTag + 1, -1);
  for (int i = 0; i < pointCount; i++) tagToIndex[nodeTags[i]] = i;

  double minX = 1e9, maxX = -1e9;
  double minY = 1e9, maxY = -1e9;
  double minZ = 1e9, maxZ = -1e9;

  for (int i = 0; i < pointCount; i++) {
    double x = nodeCoords[3 * i];
    double y = nodeCoords[3 * i + 1];
    double z = nodeCoords[3 * i + 2];

    if (x < minX) minX = x; if (x > maxX) maxX = x;
    if (y < minY) minY = y; if (y > maxY) maxY = y;
    if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
  }

  double centerX = (minX + maxX) / 2;
  double centerY = (minY + maxY) / 2;
  double centerZ = (minZ + maxZ) / 2;

  double sizeX = maxX - minX;
  double sizeY = maxY - minY;
  double sizeZ = maxZ - minZ;
  double maxSize = std::max(sizeX, std::max(sizeY, sizeZ));

  double tailStart = centerX + 0.4 * sizeX;

  int frames = 240;
  double totalTime = 3.0;
  double omega = 2 * M_PI * 1.5;
  double amplitude = 0.08 * sizeY;

  system("mkdir -p output");

  for (int f = 0; f < frames; f++) {
    double t = totalTime * f / (frames - 1);

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(pointCount);

    auto velocity = vtkSmartPointer<vtkDoubleArray>::New();
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(pointCount);

    auto wave = vtkSmartPointer<vtkDoubleArray>::New();
    wave->SetName("Wave");
    wave->SetNumberOfTuples(pointCount);

    auto color = vtkSmartPointer<vtkUnsignedCharArray>::New();
    color->SetName("Color");
    color->SetNumberOfComponents(3);
    color->SetNumberOfTuples(pointCount);

    for (int i = 0; i < pointCount; i++) {
      double x = nodeCoords[3 * i];
      double y = nodeCoords[3 * i + 1];
      double z = nodeCoords[3 * i + 2];

      double stretch = (x - centerX) / (0.5 * sizeX);
      double shiftX = amplitude * stretch * sin(omega * t);
      double shiftY = 0;

      double tail = 0;
      if (x > tailStart) tail = (x - tailStart) / (maxX - tailStart);

      shiftY += 2 * amplitude * tail * sin(omega * t);

      double vx = amplitude * stretch * omega * cos(omega * t);
      double vy = 2 * amplitude * tail * omega * cos(omega * t);

      points->SetPoint(i, x + shiftX, y + shiftY, z);
      velocity->SetTuple3(i, vx, vy, 0);

      double radius = sqrt((x - centerX) * (x - centerX) +
                           (y - centerY) * (y - centerY) +
                           (z - centerZ) * (z - centerZ));

      double scalar = sin(8 * radius / (maxSize + 1e-12) - 2 * omega * t) +
                      0.4 * cos(6 * z / (sizeZ + 1e-12) + omega * t);

      wave->SetValue(i, scalar);

      double blend = 0.5 + 0.5 * scalar;
      if (blend < 0) blend = 0; if (blend > 1) blend = 1;

      unsigned char r = (1 - blend) * 60  + blend * 220;
      unsigned char g = (1 - blend) * 0   + blend * 20;
      unsigned char b = (1 - blend) * 100 + blend * 60;

      color->SetTuple3(i, r, g, b);
    }

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);

    std::vector<std::size_t> &tetraNodes = elementNodeTags[tetraIndex];
    int tetraCount = tetraNodes.size() / 4;

    for (int i = 0; i < tetraCount; i++) {
      auto tetra = vtkSmartPointer<vtkTetra>::New();
      for (int j = 0; j < 4; j++) {
        std::size_t gmshTag = tetraNodes[4 * i + j];
        tetra->GetPointIds()->SetId(j, tagToIndex[gmshTag]);
      }
      grid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }

    grid->GetPointData()->AddArray(velocity);
    grid->GetPointData()->SetActiveVectors("Velocity");
    grid->GetPointData()->AddArray(wave);
    grid->GetPointData()->AddArray(color);
    grid->GetPointData()->SetActiveScalars("Wave");

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(("output/frame_" + std::to_string(f) + ".vtu").c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToAscii();
    writer->SetCompressor(nullptr);
    writer->Write();
  }

  gmsh::finalize();
  return 0;
}