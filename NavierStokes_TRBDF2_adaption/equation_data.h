// We start by including all the necessary deal.II header files.
//
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

// @sect{Equation data}

// In the next namespace, we declare and implement suitable functions that may be used for the initial and boundary conditions
//
namespace EquationData {
  using namespace dealii;

  static const unsigned int degree_p = 1; /*--- Polynomial degree for the pressure. The one for the velocity
                                                will be equal to degree_p + 1, but its value can be easily changed
                                                in the template parameter list of the classes with weak form ---*/

  // We declare class that describes the boundary conditions and initial one for velocity:
  //
  template<int dim>
  class Velocity: public Function<dim> {
  public:
    Velocity(const double initial_time = 0.0);

    virtual double value(const Point<dim>&  p,
                         const unsigned int component = 0) const override;

    virtual void vector_value(const Point<dim>& p,
                              Vector<double>&   values) const override;
  };


  template<int dim>
  Velocity<dim>::Velocity(const double initial_time): Function<dim>(dim, initial_time) {}


  template<int dim>
  double Velocity<dim>::value(const Point<dim>& p, const unsigned int component) const {
    AssertIndexRange(component, 3);
    if(component == 0) {
      const double Um = 1.5;
      const double H  = 4.1;

      return 4.0*Um*p(1)*(H - p(1))/(H*H);
    }
    else
      return 0.0;
  }


  template<int dim>
  void Velocity<dim>::vector_value(const Point<dim>& p, Vector<double>& values) const {
    Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));

    for(unsigned int i = 0; i < dim; ++i)
      values[i] = value(p, i);
  }


  // We do the same for the pressure
  //
  template<int dim>
  class Pressure: public Function<dim> {
  public:
    Pressure(const double initial_time = 0.0);

    virtual double value(const Point<dim>&  p,
                         const unsigned int component = 0) const override;
  };


  template<int dim>
  Pressure<dim>::Pressure(const double initial_time): Function<dim>(1, initial_time) {}


  template<int dim>
  double Pressure<dim>::value(const Point<dim>&  p, const unsigned int component) const {
    (void)component;
    AssertIndexRange(component, 1);

    return 22.0 - p(0);
  }

  // Viscosity class definition
  template<int dim, typename Number>
  class Viscosity: public Function<dim> {
  private:
    const double Cs_2 = 0.1 * 0.1;
  public:

    Viscosity(const double initial_time = 0.0);

    double value(const Point<dim>&  p, const Tensor<dim, dim, const double>& sym_grad_u, const double& dx, const double& Re) const;

    VectorizedArray<Number> value(const Point<dim, VectorizedArray<Number>>&  p, const SymmetricTensor<dim, dim, VectorizedArray<Number>>& sym_grad_u, const VectorizedArray<Number>& dx, const double& Re) const;

  };


  template<int dim, typename Number>
  Viscosity<dim, Number>::Viscosity(const double initial_time): Function<dim>(1, initial_time) {}


  template<int dim, typename Number>
  double Viscosity<dim, Number>::value(const Point<dim>&  p, const Tensor<dim, dim, const double>& sym_grad_u, const double& dx, const double& Re) const {

    // subgrid viscosity
    const double y_plus = 25.;
    const double fd = 1 - exp(-y_plus / 25.);

    double sgs_nu = Cs_2 * dx * dx * sqrt( 2 * sym_grad_u.norm());

    return 1. / Re + sgs_nu;
  }


  template<int dim, typename Number>
  VectorizedArray<Number> Viscosity<dim, Number>::value(const Point<dim, VectorizedArray<Number>>&  p, const SymmetricTensor<dim, dim, VectorizedArray<Number>>& sym_grad_u, const VectorizedArray<Number>& dx, const double& Re) const {

    // subgrid viscosity
    const double y_plus = 25.;
    const double fd = 1 - exp(-y_plus / 25.);

    auto viscosity = VectorizedArray<Number>();
    auto tensor = Tensor<dim, dim, double>();

    for(unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
      for(unsigned int di = 0; di < dim; ++di)
        for(unsigned int dj = 0; dj < dim; ++dj)
          tensor[di][dj] = sym_grad_u[di][dj][v];

      viscosity[v] = 1. / Re + Cs_2 * dx[v] * dx[v] * tensor.norm();
    }

    return viscosity;
  }

} // namespace EquationData
