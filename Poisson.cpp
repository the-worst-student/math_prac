#include <dolfin.h>
#include <mshr.h>
#include "Poisson.h"
#include <cmath>
#include <memory>

using namespace dolfin;
using namespace mshr;
class Source : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
    // Изменена правая часть: центр не в (0.5, 0.5),
    // а в верхней части круга
    double dx = x[0];
    double dy = x[1] - 0.35;
    values[0] = 8.0 * std::exp(-6.0 * (dx * dx + dy * dy));
  }
};

class dUdN : public Expression {
public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
    // Изменено граничное условие
    values[0] = 2.0 + x[0];
  }
};

class DirichletBoundary : public SubDomain {
public:
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    // Изменена область Дирихле:
    return on_boundary && x[1] >= -DOLFIN_EPS;
  }
};

int main()
{
  // Вместо квадратной области создается круг радиуса 1 с центром в начале координат
  std::shared_ptr<const CSGGeometry> domain =
      std::make_shared<Circle>(Point(0.0, 0.0), 1.0, 96);

  auto mesh = generate_mesh(domain, 48);

  auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  // Значение Дирихле изменено с 0.0 на 1.0
  auto u0 = std::make_shared<Constant>(1.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);

  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L.f = f;
  L.g = g;

  Function u(V);
  solve(a == L, u, bc);
  File file("poisson_circle.pvd");
  file << u;

  return 0;
}