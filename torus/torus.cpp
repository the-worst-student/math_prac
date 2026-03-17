#include <gmsh.h>
#include <cmath>
#include <vector>

int main() {
  gmsh::initialize();
  gmsh::model::add("torus");

  double R = 5.0;
  double r = 2.0;

  // Точки
  gmsh::model::occ::addPoint(R, 0, 0, 1.0, 1);
  gmsh::model::occ::addPoint(R, r, 0, 1.0, 2);
  gmsh::model::occ::addPoint(R, -r, 0, 1.0, 3);

  // Дуги
  gmsh::model::occ::addCircleArc(2, 1, 3, 1);
  gmsh::model::occ::addCircleArc(3, 1, 2, 2);

  // Контур и поверхность
  gmsh::model::occ::addCurveLoop({1, 2}, 1);
  gmsh::model::occ::addPlaneSurface({1}, 1);
  gmsh::model::occ::synchronize();

  // Вращение вокруг оси Y
  std::vector<std::pair<int, int>> out;
  gmsh::model::occ::revolve({{2, 1}}, 0, 0, 0, 0, 1, 0, 2 * M_PI, out);
  gmsh::model::occ::synchronize();

  // Сетка
  gmsh::model::mesh::generate(3);
  gmsh::write("torus.msh");
  gmsh::fltk::run();
  gmsh::finalize();
  return 0;
}