/* Author: Giuseppe Orlando, 2022. */

// We start by including all the necessary deal.II header files and some C++
// related ones.
//
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/fe/component_mask.h>

#include <deal.II/base/timer.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/manifold_lib.h>

#include "runtime_parameters.h"
#include "equation_data.h"

namespace MatrixFreeTools {
  using namespace dealii;

  template<int dim, typename Number, typename VectorizedArrayType>
  void compute_diagonal(const MatrixFree<dim, Number, VectorizedArrayType>&                            matrix_free,
                        LinearAlgebra::distributed::Vector<Number>&                                    diagonal_global,
                        const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType>&,
                                                 LinearAlgebra::distributed::Vector<Number>&,
                                                 const unsigned int&,
                                                 const std::pair<unsigned int, unsigned int>&)>& 	     cell_operation,
                        const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType>&,
                                                 LinearAlgebra::distributed::Vector<Number>&,
                                                 const unsigned int&,
                                                 const std::pair<unsigned int, unsigned int>&)>& 	     face_operation,
                        const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType>&,
                                                 LinearAlgebra::distributed::Vector<Number>&,
                                                 const unsigned int&,
                                                 const std::pair<unsigned int, unsigned int>&)>& 	     boundary_operation,
                        const unsigned int                                                             dof_no = 0) {
    // initialize vector
    matrix_free.initialize_dof_vector(diagonal_global, dof_no);

    const unsigned int dummy = 0;

    matrix_free.loop(cell_operation, face_operation, boundary_operation,
                     diagonal_global, dummy, false,
                     MatrixFree<dim, Number>::DataAccessOnFaces::unspecified,
                     MatrixFree<dim, Number>::DataAccessOnFaces::unspecified);
  }
}

// We include the code in a suitable namespace:
//
namespace NS_TRBDF2 {
  using namespace dealii;

  // The following class is an auxiliary one for post-processing of the vorticity
  //
  template<int dim>
  class PostprocessorVorticity: public DataPostprocessor<dim> {
  public:
    virtual void evaluate_vector_field(const DataPostprocessorInputs::Vector<dim>& inputs,
                                       std::vector<Vector<double>>&                computed_quantities) const override;

    virtual std::vector<std::string> get_names() const override;

    virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

    virtual UpdateFlags get_needed_update_flags() const override;
  };

  // This function evaluates the vorticty in both 2D and 3D cases
  //
  template <int dim>
  void PostprocessorVorticity<dim>::evaluate_vector_field(const DataPostprocessorInputs::Vector<dim>& inputs,
                                                          std::vector<Vector<double>>&                computed_quantities) const {
    const unsigned int n_quadrature_points = inputs.solution_values.size();

    /*--- Check the correctness of all data structres ---*/
    Assert(inputs.solution_gradients.size() == n_quadrature_points, ExcInternalError());
    Assert(computed_quantities.size() == n_quadrature_points, ExcInternalError());

    Assert(inputs.solution_values[0].size() == dim, ExcInternalError());

    if(dim == 2) {
      Assert(computed_quantities[0].size() == 1, ExcInternalError());
    }
    else {
      Assert(computed_quantities[0].size() == dim, ExcInternalError());
    }

    /*--- Compute the vorticty ---*/
    if(dim == 2) {
      for(unsigned int q = 0; q < n_quadrature_points; ++q)
        computed_quantities[q](0) = inputs.solution_gradients[q][1][0] - inputs.solution_gradients[q][0][1];
    }
    else {
      for(unsigned int q = 0; q < n_quadrature_points; ++q) {
        computed_quantities[q](0) = inputs.solution_gradients[q][2][1] - inputs.solution_gradients[q][1][2];
        computed_quantities[q](1) = inputs.solution_gradients[q][0][2] - inputs.solution_gradients[q][2][0];
        computed_quantities[q](2) = inputs.solution_gradients[q][1][0] - inputs.solution_gradients[q][0][1];
      }
    }
  }

  // This auxiliary function is required by the base class DataProcessor and simply
  // sets the name for the output file
  //
  template<int dim>
  std::vector<std::string> PostprocessorVorticity<dim>::get_names() const {
    std::vector<std::string> names;
    names.emplace_back("vorticity");
    if(dim == 3) {
      names.emplace_back("vorticity");
      names.emplace_back("vorticity");
    }

    return names;
  }

  // This auxiliary function is required by the base class DataProcessor and simply
  // specifies if the vorticity is a scalar (2D) or a vector (3D)
  //
  template<int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  PostprocessorVorticity<dim>::get_data_component_interpretation() const {
    std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
    if(dim == 2)
      interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    else {
      interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
      interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
      interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    }

    return interpretation;
  }

  // This auxiliary function is required by the base class DataProcessor and simply
  // sets which variables have to updated (only the gradients)
  //
  template<int dim>
  UpdateFlags PostprocessorVorticity<dim>::get_needed_update_flags() const {
    return update_gradients;
  }


  // The following structs are auxiliary objects for mesh refinement. ScratchData simply sets
  // the FEValues object
  //
  template <int dim>
  struct ScratchData {
    ScratchData(const FiniteElement<dim>& fe,
                const unsigned int        quadrature_degree,
                const UpdateFlags         update_flags): fe_values(fe, QGauss<dim>(quadrature_degree), update_flags) {}

    ScratchData(const ScratchData<dim>& scratch_data): fe_values(scratch_data.fe_values.get_fe(),
                                                                 scratch_data.fe_values.get_quadrature(),
                                                                 scratch_data.fe_values.get_update_flags()) {}
    FEValues<dim> fe_values;
  };


  // CopyData simply sets the cell index
  //
  struct CopyData {
    CopyData() : cell_index(numbers::invalid_unsigned_int), value(0.0) {}

    CopyData(const CopyData &) = default;

    unsigned int cell_index;
    double       value;
  };


  // @sect{ <code>NavierStokesProjectionOperator::NavierStokesProjectionOperator</code> }

  // The following class sets effecively the weak formulation of the problems for the different stages
  // and for both velocity and pressure.
  // The template parameters are the dimnesion of the problem, the polynomial degree for the pressure,
  // the polynomial degree for the velocity, the number of quadrature points for integrals for the pressure step,
  // the number of quadrature points for integrals for the velocity step, the type of vector for storage and the type
  // of floating point data (in general double or float for preconditioners structures if desired).
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  class NavierStokesProjectionOperator: public MatrixFreeOperators::Base<dim, Vec> {
  public:
    using Number = typename Vec::value_type;

    NavierStokesProjectionOperator();

    NavierStokesProjectionOperator(RunTimeParameters::Data_Storage& data);

    void set_dt(const double time_step);

    void set_Reynolds(const double reynolds);

    void set_TR_BDF2_stage(const unsigned int stage);

    void set_NS_stage(const unsigned int stage);

    void set_u_extr(const Vec& src);

    void set_deltas(const Vec& src);

    void set_y_plus(const Vec& src);

    void vmult_rhs_velocity(Vec& dst, const std::vector<Vec>& src) const;

    void vmult_rhs_pressure(Vec& dst, const std::vector<Vec>& src) const;

    void vmult_grad_p_projection(Vec& dst, const Vec& src) const;

    virtual void compute_diagonal() override;

  protected:
    double       Re;
    double       dt;
    bool         no_slip;

    /*--- Parameters of time-marching scheme ---*/
    double       gamma;
    double       a31;
    double       a32;
    double       a33;

    unsigned int TR_BDF2_stage; /*--- Flag to denote at which stage of the TR-BDF2 are ---*/
    unsigned int NS_stage;      /*--- Flag to denote at which stage of NS solution inside each TR-BDF2 stage we are
                                      (solution of the velocity or of the pressure)---*/

    virtual void apply_add(Vec& dst, const Vec& src) const override;

  private:
    /*--- Auxiliary variable for the TR stage
          (just to avoid to report a lot of 0.5 and for my personal choice to be coherent with the article) ---*/
    const double a21 = 0.5;
    const double a22 = 0.5;

    /*--- Penalty method parameters, theta = 1 means SIP, while C_p and C_u are the penalization coefficients ---*/
    const double theta_v = 1.0;
    const double theta_p = 1.0;
    const double C_p = 1.0*(fe_degree_p + 1)*(fe_degree_p + 1);
    const double C_u = 1.0*(fe_degree_v + 1)*(fe_degree_v + 1);

    Vec u_extr; /*--- Auxiliary variable to update the extrapolated velocity ---*/
    Vec deltas;
    Vec y_plus;

    EquationData::Velocity<dim>  vel_boundary_inflow; /*--- Auxiliary variable to impose velocity boundary conditions ---*/

    EquationData::Viscosity<dim, Number> viscosity;

    /*--- The following functions basically assemble the linear and bilinear forms. Their syntax is due to
          the base class MatrixFreeOperators::Base ---*/
    void assemble_rhs_cell_term_velocity(const MatrixFree<dim, Number>&               data,
                                         Vec&                                         dst,
                                         const std::vector<Vec>&                      src,
                                         const std::pair<unsigned int, unsigned int>& cell_range) const;
    void assemble_rhs_face_term_velocity(const MatrixFree<dim, Number>&               data,
                                         Vec&                                         dst,
                                         const std::vector<Vec>&                      src,
                                         const std::pair<unsigned int, unsigned int>& face_range) const;
    void assemble_rhs_boundary_term_velocity(const MatrixFree<dim, Number>&               data,
                                             Vec&                                         dst,
                                             const std::vector<Vec>&                      src,
                                             const std::pair<unsigned int, unsigned int>& face_range) const;

    void assemble_rhs_cell_term_pressure(const MatrixFree<dim, Number>&               data,
                                         Vec&                                         dst,
                                         const std::vector<Vec>&                      src,
                                         const std::pair<unsigned int, unsigned int>& cell_range) const;
    void assemble_rhs_face_term_pressure(const MatrixFree<dim, Number>&               data,
                                         Vec&                                         dst,
                                         const std::vector<Vec>&                      src,
                                         const std::pair<unsigned int, unsigned int>& face_range) const;
    void assemble_rhs_boundary_term_pressure(const MatrixFree<dim, Number>&               data,
                                             Vec&                                         dst,
                                             const std::vector<Vec>&                      src,
                                             const std::pair<unsigned int, unsigned int>& face_range) const;

    void assemble_cell_term_velocity(const MatrixFree<dim, Number>&               data,
                                     Vec&                                         dst,
                                     const Vec&                                   src,
                                     const std::pair<unsigned int, unsigned int>& cell_range) const;
    void assemble_face_term_velocity(const MatrixFree<dim, Number>&               data,
                                     Vec&                                         dst,
                                     const Vec&                                   src,
                                     const std::pair<unsigned int, unsigned int>& face_range) const;
    void assemble_boundary_term_velocity(const MatrixFree<dim, Number>&               data,
                                         Vec&                                         dst,
                                         const Vec&                                   src,
                                         const std::pair<unsigned int, unsigned int>& face_range) const;

    void assemble_cell_term_pressure(const MatrixFree<dim, Number>&               data,
                                     Vec&                                         dst,
                                     const Vec&                                   src,
                                     const std::pair<unsigned int, unsigned int>& cell_range) const;
    void assemble_face_term_pressure(const MatrixFree<dim, Number>&               data,
                                     Vec&                                         dst,
                                     const Vec&                                   src,
                                     const std::pair<unsigned int, unsigned int>& face_range) const;
    void assemble_boundary_term_pressure(const MatrixFree<dim, Number>&               data,
                                         Vec&                                         dst,
                                         const Vec&                                   src,
                                         const std::pair<unsigned int, unsigned int>& face_range) const;

    void assemble_cell_term_projection_grad_p(const MatrixFree<dim, Number>&               data,
                                              Vec&                                         dst,
                                              const Vec&                                   src,
                                              const std::pair<unsigned int, unsigned int>& cell_range) const;
    void assemble_rhs_cell_term_projection_grad_p(const MatrixFree<dim, Number>&               data,
                                                  Vec&                                         dst,
                                                  const Vec&                                   src,
                                                  const std::pair<unsigned int, unsigned int>& cell_range) const;

    void assemble_diagonal_cell_term_velocity(const MatrixFree<dim, Number>&               data,
                                              Vec&                                         dst,
                                              const unsigned int&                          src,
                                              const std::pair<unsigned int, unsigned int>& cell_range) const;
    void assemble_diagonal_face_term_velocity(const MatrixFree<dim, Number>&               data,
                                              Vec&                                         dst,
                                              const unsigned int&                          src,
                                              const std::pair<unsigned int, unsigned int>& face_range) const;
    void assemble_diagonal_boundary_term_velocity(const MatrixFree<dim, Number>&               data,
                                                  Vec&                                         dst,
                                                  const unsigned int&                          src,
                                                  const std::pair<unsigned int, unsigned int>& face_range) const;

    void assemble_diagonal_cell_term_pressure(const MatrixFree<dim, Number>&               data,
                                              Vec&                                         dst,
                                              const unsigned int&                          src,
                                              const std::pair<unsigned int, unsigned int>& cell_range) const;
    void assemble_diagonal_face_term_pressure(const MatrixFree<dim, Number>&               data,
                                              Vec&                                         dst,
                                              const unsigned int&                          src,
                                              const std::pair<unsigned int, unsigned int>& face_range) const;
    void assemble_diagonal_boundary_term_pressure(const MatrixFree<dim, Number>&               data,
                                                  Vec&                                         dst,
                                                  const unsigned int&                          src,
                                                  const std::pair<unsigned int, unsigned int>& face_range) const;
  };


  // We start with the default constructor. It is important for MultiGrid, so it is fundamental
  // to properly set the parameters of the time scheme.
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  NavierStokesProjectionOperator():
    MatrixFreeOperators::Base<dim, Vec>(), Re(), dt(), no_slip(true), gamma(2.0 - std::sqrt(2.0)), a31((1.0 - gamma)/(2.0*(2.0 - gamma))),
                                           a32(a31), a33(1.0/(2.0 - gamma)), TR_BDF2_stage(1), NS_stage(1), u_extr(), deltas(), y_plus() {}


  // We focus now on the constructor with runtime parameters storage
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  NavierStokesProjectionOperator(RunTimeParameters::Data_Storage& data):
    MatrixFreeOperators::Base<dim, Vec>(), Re(data.Reynolds), dt(data.dt), no_slip(data.no_slip),
                                           gamma(2.0 - std::sqrt(2.0)), a31((1.0 - gamma)/(2.0*(2.0 - gamma))),
                                           a32(a31), a33(1.0/(2.0 - gamma)), TR_BDF2_stage(1), NS_stage(1), u_extr(), deltas(), y_plus(),
                                           vel_boundary_inflow(data.initial_time),
                                           viscosity(data.initial_time, data.Cs2) {}


  // Setter of time-step (called by Multigrid and in case a smaller time-step towards the end is needed)
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  set_dt(const double time_step) {
    dt = time_step;
  }


  // Setter of Reynolds number
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  set_Reynolds(const double reynolds) {
    Re = reynolds;
  }


  // Setter of TR-BDF2 stage (this can be known only during the effective execution
  // and so it has to be demanded to the class that really solves the problem)
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  set_TR_BDF2_stage(const unsigned int stage) {
    AssertIndexRange(stage, 3);
    Assert(stage > 0, ExcInternalError());

    TR_BDF2_stage = stage;
  }


  // Setter of NS stage (this can be known only during the effective execution
  // and so it has to be demanded to the class that really solves the problem)
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  set_NS_stage(const unsigned int stage) {
    AssertIndexRange(stage, 4);
    Assert(stage > 0, ExcInternalError());

    NS_stage = stage;
  }


  // Setter of extrapolated velocity for different stages
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  set_u_extr(const Vec& src) {
    u_extr = src;
    u_extr.update_ghost_values();
  }

  // Setter of deltas
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  set_deltas(const Vec& src) {
    deltas = src;
    deltas.update_ghost_values();
  }

  // Setter of y_plus
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  set_y_plus(const Vec& src) {
    y_plus = src;
    y_plus.update_ghost_values();
  }


  // We are in a DG-MatrixFree framework, so it is convenient to compute separately cell contribution,
  // internal faces contributions and boundary faces contributions. We start by
  // assembling the rhs cell term for the velocity.
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_rhs_cell_term_velocity(const MatrixFree<dim, Number>&               data,
                                  Vec&                                         dst,
                                  const std::vector<Vec>&                      src,
                                  const std::pair<unsigned int, unsigned int>& cell_range) const {
    if(TR_BDF2_stage == 1) {
      /*--- We first start by declaring the suitable instances to read the old velocity, the
      extrapolated velocity and the old pressure. 'phi' will be used only to submit the result.
      The second argument specifies which dof handler has to be used (in this implementation 0 stands for
      velocity and 1 for pressure). ---*/
      FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0),
                                                                   phi_old(data, 0),
                                                                   phi_old_extr(data, 0),
                                                                   phi_force(data, 0);
      FEEvaluation<dim, fe_degree_p, n_q_points_1d_v, 1, Number>   phi_old_press(data, 1);
      FEEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, 2),
                                                                   phi_y_plus(data, 2);

      /*--- We loop over the cells in the range ---*/
      for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
        /*--- Now we need to assign the current cell to each FEEvaluation object and then to specify which src vector
        it has to read (the proper order is clearly delegated to the user, which has to pay attention in the function
        call to be coherent). ---*/
        phi_old.reinit(cell);
        phi_old.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
                                                           /*--- The 'gather_evaluate' function reads data from the vector.
                                                           The second and third parameter specifies if you want to read
                                                           values and/or derivative related quantities ---*/
        phi_old_extr.reinit(cell);
        phi_old_extr.gather_evaluate(src[1], EvaluationFlags::values);
        phi_old_press.reinit(cell);
        phi_old_press.gather_evaluate(src[2], EvaluationFlags::values);
        phi.reinit(cell);

        phi_deltas.reinit(cell);
        phi_deltas.gather_evaluate(src[3], EvaluationFlags::values);

        phi_force.reinit(cell);
        phi_force.gather_evaluate(src[4], EvaluationFlags::values);

        phi_y_plus.reinit(cell);
        phi_y_plus.gather_evaluate(src[5], EvaluationFlags::values);

        /*--- Now we loop over all the quadrature points to compute the integrals ---*/
        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          const auto& u_n                = phi_old.get_value(q);
          const auto& grad_u_n           = 2.0 * phi_old.get_symmetric_gradient(q);

          const auto& u_n_gamma_ov_2     = phi_old_extr.get_value(q);
          const auto& tensor_product_u_n = outer_product(u_n, u_n_gamma_ov_2);
          const auto& p_n                = phi_old_press.get_value(q);

          auto p_n_times_identity        = tensor_product_u_n;
          p_n_times_identity = 0;
          for(unsigned int d = 0; d < dim; ++d)
            p_n_times_identity[d][d] = p_n;

          const auto& dx                 = phi_deltas.get_value(q);

          phi.submit_value(1.0/(gamma*dt)*u_n, q); 

          // phi.submit_value(1.0/(gamma*dt)*u_n + phi_force.get_value(q), q); /*--- 'submit_value' contains quantites that we want to test against the
          //                                                 test function ---*/
          phi.submit_gradient(-a21*viscosity.value(phi_y_plus.get_value(q), grad_u_n, dx, Re)*grad_u_n +
                               a21*tensor_product_u_n + p_n_times_identity, q);
          /*--- 'submit_gradient' contains quantites that we want to test against the gradient of test function ---*/
        }
        phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
        /*--- 'integrate_scatter' is the responsible of distributing into dst.
              The flag parameter specifies if we are testing against the test function and/or its gradient ---*/
      }
    }
    else {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0),
                                                                   phi_old(data, 0),
                                                                   phi_int(data, 0),
                                                                   phi_force(data, 0);
      FEEvaluation<dim, fe_degree_p, n_q_points_1d_v, 1, Number>   phi_old_press(data, 1);
      FEEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, 2),
                                                                   phi_y_plus(data, 2);

      /*--- We loop over the cells in the range ---*/
      for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
        phi_old.reinit(cell);
        phi_old.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_int.reinit(cell);
        phi_int.gather_evaluate(src[1], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_press.reinit(cell);
        phi_old_press.gather_evaluate(src[2], EvaluationFlags::values);
        phi.reinit(cell);

        phi_deltas.reinit(cell);
        phi_deltas.gather_evaluate(src[4], EvaluationFlags::values);

        phi_force.reinit(cell);
        phi_force.gather_evaluate(src[5], EvaluationFlags::values);

        phi_y_plus.reinit(cell);
        phi_y_plus.gather_evaluate(src[6], EvaluationFlags::values);

        /*--- Now we loop over all the quadrature points to compute the integrals ---*/
        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          const auto& u_n                      = phi_old.get_value(q);
          const auto& grad_u_n                 = 2.0 * phi_old.get_symmetric_gradient(q);

          const auto& u_n_gamma                = phi_int.get_value(q);
          const auto& grad_u_n_gamma           = 2.0 * phi_int.get_symmetric_gradient(q);
          const auto& tensor_product_u_n       = outer_product(u_n, u_n);
          const auto& tensor_product_u_n_gamma = outer_product(u_n_gamma, u_n_gamma);

          const auto& p_n                      = phi_old_press.get_value(q);
          auto p_n_times_identity              = tensor_product_u_n;
          p_n_times_identity = 0;
          for(unsigned int d = 0; d < dim; ++d)
            p_n_times_identity[d][d] = p_n;

          const auto& dx                       = phi_deltas.get_value(q);

          phi.submit_value(1.0/((1.0 - gamma)*dt)*u_n_gamma + phi_force.get_value(q), q);
          phi.submit_gradient(a32*tensor_product_u_n_gamma + a31*tensor_product_u_n -
                             a32*viscosity.value(phi_y_plus.get_value(q), grad_u_n_gamma, dx, Re)*grad_u_n_gamma -
                             a31*viscosity.value(phi_y_plus.get_value(q), grad_u_n_gamma, dx, Re)*grad_u_n + p_n_times_identity, q);

        }
        phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
  }


  // The followinf function assembles rhs face term for the velocity
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_rhs_face_term_velocity(const MatrixFree<dim, Number>&               data,
                                  Vec&                                         dst,
                                  const std::vector<Vec>&                      src,
                                  const std::pair<unsigned int, unsigned int>& face_range) const {
    if(TR_BDF2_stage == 1) {
      /*--- We first start by declaring the suitable instances to read already available quantities. In this case
      we are at the face between two elements and this is the reason of 'FEFaceEvaluation'. It contains an extra
      input argument, the second one, that specifies if it is from 'interior' or not---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi_p(data, true, 0),
                                                                       phi_m(data, false, 0),
                                                                       phi_old_p(data, true, 0),
                                                                       phi_old_m(data, false, 0),
                                                                       phi_old_extr_p(data, true, 0),
                                                                       phi_old_extr_m(data, false, 0);
      FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_v, 1, Number>   phi_old_press_p(data, true, 1),
                                                                       phi_old_press_m(data, false, 1);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas_p(data, true, 2),
                                                                       phi_deltas_m(data, false, 2),
                                                                       phi_y_plus_p(data, true, 2),
                                                                       phi_y_plus_m(data, false, 2);

      /*--- We loop over the faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_old_p.reinit(face);
        phi_old_p.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_m.reinit(face);
        phi_old_m.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_extr_p.reinit(face);
        phi_old_extr_p.gather_evaluate(src[1], EvaluationFlags::values);
        phi_old_extr_m.reinit(face);
        phi_old_extr_m.gather_evaluate(src[1], EvaluationFlags::values);
        phi_old_press_p.reinit(face);
        phi_old_press_p.gather_evaluate(src[2], EvaluationFlags::values);
        phi_old_press_m.reinit(face);
        phi_old_press_m.gather_evaluate(src[2], EvaluationFlags::values);
        phi_p.reinit(face);
        phi_m.reinit(face);

        phi_deltas_p.reinit(face);
        phi_deltas_p.gather_evaluate(src[3], EvaluationFlags::values);
        phi_deltas_m.reinit(face);
        phi_deltas_m.gather_evaluate(src[3], EvaluationFlags::values);

        phi_y_plus_p.reinit(face);
        phi_y_plus_p.gather_evaluate(src[5], EvaluationFlags::values);
        phi_y_plus_m.reinit(face);
        phi_y_plus_m.gather_evaluate(src[5], EvaluationFlags::values);

        /*--- Now we loop over all the quadrature points to compute the integrals ---*/
        for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
          const auto& n_plus                 = phi_p.get_normal_vector(q); /*--- The normal vector is the same
                                                                                 for both phi_p and phi_m. If the face is interior,
                                                                                 it correspond to the outer normal ---*/
          const auto& dx_p                   = phi_deltas_p.get_value(q);
          const auto& dx_m                   = phi_deltas_m.get_value(q);

          const auto& avg_visc_grad_u_old    = viscosity.value(phi_y_plus_p.get_value(q), phi_old_p.get_symmetric_gradient(q), dx_p, Re)*
                                               phi_old_p.get_symmetric_gradient(q) +
                                               viscosity.value(phi_y_plus_m.get_value(q), phi_old_m.get_symmetric_gradient(q), dx_m, Re)*
                                               phi_old_m.get_symmetric_gradient(q);
          const auto& avg_tensor_product_u_n = 0.5*(outer_product(phi_old_p.get_value(q), phi_old_extr_p.get_value(q)) +
                                                    outer_product(phi_old_m.get_value(q), phi_old_extr_m.get_value(q)));
          const auto& avg_p_old              = 0.5*(phi_old_press_p.get_value(q) + phi_old_press_m.get_value(q));

          phi_p.submit_value((a21*avg_visc_grad_u_old - a21*avg_tensor_product_u_n)*n_plus - avg_p_old*n_plus, q);
          phi_m.submit_value(-(a21*avg_visc_grad_u_old - a21*avg_tensor_product_u_n)*n_plus + avg_p_old*n_plus, q);

        }
        phi_p.integrate_scatter(EvaluationFlags::values, dst);
        phi_m.integrate_scatter(EvaluationFlags::values, dst);
      }
    }
    else {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi_p(data, true, 0),
                                                                       phi_m(data, false, 0),
                                                                       phi_old_p(data, true, 0),
                                                                       phi_old_m(data, false, 0),
                                                                       phi_int_p(data, true, 0),
                                                                       phi_int_m(data, false, 0);
      FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_v, 1, Number>   phi_old_press_p(data, true, 1),
                                                                       phi_old_press_m(data, false, 1);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas_p(data, true, 2),
                                                                       phi_deltas_m(data, false, 2),
                                                                       phi_y_plus_p(data, true, 2),
                                                                       phi_y_plus_m(data, false, 2);

      /*--- We loop over the faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++ face) {
        phi_old_p.reinit(face);
        phi_old_p.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_m.reinit(face);
        phi_old_m.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_int_p.reinit(face);
        phi_int_p.gather_evaluate(src[1], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_int_m.reinit(face);
        phi_int_m.gather_evaluate(src[1], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_press_p.reinit(face);
        phi_old_press_p.gather_evaluate(src[2], EvaluationFlags::values);
        phi_old_press_m.reinit(face);
        phi_old_press_m.gather_evaluate(src[2], EvaluationFlags::values);
        phi_p.reinit(face);
        phi_m.reinit(face);

        phi_deltas_p.reinit(face);
        phi_deltas_p.gather_evaluate(src[4], EvaluationFlags::values);
        phi_deltas_m.reinit(face);
        phi_deltas_m.gather_evaluate(src[4], EvaluationFlags::values);

        phi_y_plus_p.reinit(face);
        phi_y_plus_p.gather_evaluate(src[6], EvaluationFlags::values);
        phi_y_plus_m.reinit(face);
        phi_y_plus_m.gather_evaluate(src[6], EvaluationFlags::values);

        /*--- Now we loop over all the quadrature points to compute the integrals ---*/
        for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
          const auto& n_plus                      = phi_p.get_normal_vector(q);

          const auto& dx_p                        = phi_deltas_p.get_value(q);
          const auto& dx_m                        = phi_deltas_m.get_value(q);

          const auto& avg_visc_grad_u_old         = viscosity.value(phi_y_plus_p.get_value(q), phi_old_p.get_symmetric_gradient(q), dx_p, Re)*
                                                    phi_old_p.get_symmetric_gradient(q) +
                                                    viscosity.value(phi_y_plus_m.get_value(q), phi_old_m.get_symmetric_gradient(q), dx_m, Re)*
                                                    phi_old_m.get_symmetric_gradient(q);

          const auto& avg_visc_grad_u_int         = viscosity.value(phi_y_plus_p.get_value(q), phi_int_p.get_symmetric_gradient(q), dx_p, Re)*
                                                    phi_int_p.get_symmetric_gradient(q) +
                                                    viscosity.value(phi_y_plus_m.get_value(q), phi_int_m.get_symmetric_gradient(q), dx_m, Re)*
                                                    phi_int_m.get_symmetric_gradient(q);

          const auto& avg_tensor_product_u_n       = 0.5*(outer_product(phi_old_p.get_value(q), phi_old_p.get_value(q)) +
                                                          outer_product(phi_old_m.get_value(q), phi_old_m.get_value(q)));
          const auto& avg_tensor_product_u_n_gamma = 0.5*(outer_product(phi_int_p.get_value(q), phi_int_p.get_value(q)) +
                                                          outer_product(phi_int_m.get_value(q), phi_int_m.get_value(q)));
          const auto& avg_p_old                    = 0.5*(phi_old_press_p.get_value(q) + phi_old_press_m.get_value(q));


          phi_p.submit_value((a31*avg_visc_grad_u_old + a32*avg_visc_grad_u_int -
                              a31*avg_tensor_product_u_n - a32*avg_tensor_product_u_n_gamma)*n_plus - avg_p_old*n_plus, q);
          phi_m.submit_value(-(a31*avg_visc_grad_u_old + a32*avg_visc_grad_u_int -
                              a31*avg_tensor_product_u_n - a32*avg_tensor_product_u_n_gamma)*n_plus + avg_p_old*n_plus, q);
        }
        phi_p.integrate_scatter(EvaluationFlags::values, dst);
        phi_m.integrate_scatter(EvaluationFlags::values, dst);
      }
    }
  }


  // The followinf function assembles rhs boundary term for the velocity
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_rhs_boundary_term_velocity(const MatrixFree<dim, Number>&               data,
                                      Vec&                                         dst,
                                      const std::vector<Vec>&                      src,
                                      const std::pair<unsigned int, unsigned int>& face_range) const {
    if(TR_BDF2_stage == 1) {
      /*--- We first start by declaring the suitable instances to read already available quantities. Clearly on the boundary
      the second argument has to be true. ---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, true, 0),
                                                                       phi_old(data, true, 0),
                                                                       phi_old_extr(data, true, 0);
      FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_v, 1, Number>   phi_old_press(data, true, 1);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, true, 2),
                                                                       phi_y_plus(data, true, 2);

      /*--- We loop over the faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_old.reinit(face);
        phi_old.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_extr.reinit(face);
        phi_old_extr.gather_evaluate(src[1], EvaluationFlags::values);
        phi_old_press.reinit(face);
        phi_old_press.gather_evaluate(src[2], EvaluationFlags::values);
        phi.reinit(face);

        phi_deltas.reinit(face);
        phi_deltas.gather_evaluate(src[3], EvaluationFlags::values);

        phi_y_plus.reinit(face);
        phi_y_plus.gather_evaluate(src[5], EvaluationFlags::values);

        const auto boundary_id = data.get_boundary_id(face); /*--- Get the id in order to impose the proper boundary condition ---*/
        
        const auto coef_jump   = (boundary_id == 1 || (!no_slip && boundary_id == 3)) ? 0.0 : C_u*std::abs((phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]);
        const double aux_coeff = (boundary_id == 1 || (!no_slip && boundary_id == 3)) ? 0.0 : 1.0;

        /*--- Now we loop over all the quadrature points to compute the integrals ---*/
        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          const auto& n_plus             = phi.get_normal_vector(q);

          const auto& grad_u_old         = 2.0 * phi_old.get_symmetric_gradient(q);
          const auto& tensor_product_u_n = outer_product(phi_old.get_value(q), phi_old_extr.get_value(q));
          const auto& p_old              = phi_old_press.get_value(q);

          const auto& dx                 = phi_deltas.get_value(q);

          const auto& point_vectorized   = phi.quadrature_point(q);
          auto u_int_m                   = Tensor<1, dim, VectorizedArray<Number>>();
          if(boundary_id == 0) {
            for(unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
              Point<dim> point; /*--- The point returned by the 'quadrature_point' function is not an instance of Point
                                      and so it is not ready to be directly used. We need to pay attention to the
                                      vectorization ---*/
              for(unsigned int d = 0; d < dim; ++d){
                point[d] = point_vectorized[d][v];
              }
              for(unsigned int d = 0; d < dim; ++d)
                u_int_m[d][v] = vel_boundary_inflow.value(point, d);

            }
          }
          const auto& tensor_product_u_int_m = outer_product(u_int_m, phi_old_extr.get_value(q));
          const auto& lambda                 = (boundary_id == 1 || (!no_slip && boundary_id == 3)) ?
                                               0.0 : std::abs(scalar_product(phi_old_extr.get_value(q), n_plus));                          

          phi.submit_value((a21*viscosity.value(phi_y_plus.get_value(q), grad_u_old, dx, Re)*grad_u_old - a21*tensor_product_u_n)*n_plus -
                           p_old*n_plus + a22*2.0*viscosity.value(phi_y_plus.get_value(q), grad_u_old, dx, Re)*coef_jump*u_int_m -
                           aux_coeff*a22*tensor_product_u_int_m*n_plus + a22*lambda*u_int_m, q);
          phi.submit_gradient(-aux_coeff*theta_v*a22*viscosity.value(phi_y_plus.get_value(q), grad_u_old, dx, Re)*
                             (outer_product(u_int_m, n_plus) + outer_product(n_plus, u_int_m)), q);
        }
        phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
    else {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, true, 0),
                                                                       phi_old(data, true, 0),
                                                                       phi_int(data, true, 0),
                                                                       phi_int_extr(data, true, 0);
      FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_v, 1, Number>   phi_old_press(data, true, 1);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, true, 2),
                                                                       phi_y_plus(data, true, 2);

      /*--- We loop over the faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++ face) {
        phi_old.reinit(face);
        phi_old.gather_evaluate(src[0], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_int.reinit(face);
        phi_int.gather_evaluate(src[1], EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_press.reinit(face);
        phi_old_press.gather_evaluate(src[2], EvaluationFlags::values);
        phi_int_extr.reinit(face);
        phi_int_extr.gather_evaluate(src[3], EvaluationFlags::values);
        phi.reinit(face);

        phi_deltas.reinit(face);
        phi_deltas.gather_evaluate(src[4], EvaluationFlags::values);

        phi_y_plus.reinit(face);
        phi_y_plus.gather_evaluate(src[6], EvaluationFlags::values);

        const auto boundary_id = data.get_boundary_id(face);
        const auto coef_jump   = (boundary_id == 1 || (!no_slip && boundary_id == 3)) ?
                                 0.0 : C_u*std::abs((phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]);
        const double aux_coeff = (boundary_id == 1 || (!no_slip && boundary_id == 3)) ? 0.0 : 1.0;

        /*--- Now we loop over all the quadrature points to compute the integrals ---*/
        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          const auto& n_plus                   = phi.get_normal_vector(q);

          const auto& grad_u_old               = 2. * phi_old.get_symmetric_gradient(q);
          const auto& grad_u_int               = 2. * phi_int.get_symmetric_gradient(q);
          const auto& tensor_product_u_n       = outer_product(phi_old.get_value(q), phi_old.get_value(q));
          const auto& tensor_product_u_n_gamma = outer_product(phi_int.get_value(q), phi_int.get_value(q));
          const auto& p_old                    = phi_old_press.get_value(q);
          const VectorizedArray<Number>& dx = phi_deltas.get_value(q);
          const auto& point_vectorized         = phi.quadrature_point(q);
          auto u_m                             = Tensor<1, dim, VectorizedArray<Number>>();
          if(boundary_id == 0) {
            for(unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
              Point<dim> point;
              for(unsigned int d = 0; d < dim; ++d)
                point[d] = point_vectorized[d][v];
              for(unsigned int d = 0; d < dim; ++d)
                u_m[d][v] = vel_boundary_inflow.value(point, d);
            }
          }
          const auto& tensor_product_u_m = outer_product(u_m, phi_int_extr.get_value(q));
          const auto& lambda             = (boundary_id == 1 || (!no_slip && boundary_id == 3)) ?
                                           0.0 : std::abs(scalar_product(phi_int_extr.get_value(q), n_plus));

          const auto& visc_int           = viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re);
          const auto& visc_old           = viscosity.value(phi_y_plus.get_value(q), grad_u_old, dx, Re);

          phi.submit_value((a31*visc_old*grad_u_old + a32*visc_int*grad_u_int -
                          a31*tensor_product_u_n - a32*tensor_product_u_n_gamma)*n_plus - p_old*n_plus +
                          a33*visc_int*2.0*coef_jump*u_m -
                          aux_coeff*a33*tensor_product_u_m*n_plus + a33*lambda*u_m, q);
          phi.submit_gradient(-aux_coeff*theta_v*a33*visc_int*(outer_product(u_m, n_plus) + outer_product(n_plus, u_m)), q);
        }
        phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
  }


  // Put together all the previous steps for velocity. This is done automatically by the loop function of 'MatrixFree' class
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  vmult_rhs_velocity(Vec& dst, const std::vector<Vec>& src) const {
    for(auto& vec : src)
      vec.update_ghost_values();

    this->data->loop(&NavierStokesProjectionOperator::assemble_rhs_cell_term_velocity,
                     &NavierStokesProjectionOperator::assemble_rhs_face_term_velocity,
                     &NavierStokesProjectionOperator::assemble_rhs_boundary_term_velocity,
                     this, dst, src, true,
                     MatrixFree<dim, Number>::DataAccessOnFaces::unspecified,
                     MatrixFree<dim, Number>::DataAccessOnFaces::unspecified);
  }


  // Now we focus on computing the rhs for the projection step for the pressure with the same ratio.
  // The following function assembles rhs cell term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_rhs_cell_term_pressure(const MatrixFree<dim, Number>&               data,
                                  Vec&                                         dst,
                                  const std::vector<Vec>&                      src,
                                  const std::pair<unsigned int, unsigned int>& cell_range) const {
    /*--- We first start by declaring the suitable instances to read already available quantities.
          The third parameter specifies that we want to use the second quadrature formula stored. ---*/
    FEEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number>   phi(data, 1, 1),
                                                                 phi_old(data, 1, 1);
    FEEvaluation<dim, fe_degree_v, n_q_points_1d_p, dim, Number> phi_proj(data, 0, 1);

    const double coeff   = (TR_BDF2_stage == 1) ? 1.0e6*gamma*dt*gamma*dt : 1.0e6*(1.0 - gamma)*dt*(1.0 - gamma)*dt;

    const double coeff_2 = (TR_BDF2_stage == 1) ? gamma*dt : (1.0 - gamma)*dt;

    /*--- We loop over cells in the range ---*/
    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
      phi_proj.reinit(cell);
      phi_proj.gather_evaluate(src[0], EvaluationFlags::values);
      phi_old.reinit(cell);
      phi_old.gather_evaluate(src[1], EvaluationFlags::values);
      phi.reinit(cell);

      /*--- Now we loop over all the quadrature points to compute the integrals ---*/
      for(unsigned int q = 0; q < phi.n_q_points; ++q) {
        const auto& u_star_star = phi_proj.get_value(q);
        const auto& p_old       = phi_old.get_value(q);

        phi.submit_value(1.0/coeff*p_old, q);
        phi.submit_gradient(1.0/coeff_2*u_star_star, q);
      }
      phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
    }
  }


  // The following function assembles rhs face term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_rhs_face_term_pressure(const MatrixFree<dim, Number>&               data,
                                  Vec&                                         dst,
                                  const std::vector<Vec>&                      src,
                                  const std::pair<unsigned int, unsigned int>& face_range) const {
    /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
    FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number>   phi_p(data, true, 1, 1),
                                                                     phi_m(data, false, 1, 1);
    FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_p, dim, Number> phi_proj_p(data, true, 0, 1),
                                                                     phi_proj_m(data, false, 0, 1);

    const double coeff = (TR_BDF2_stage == 1) ? 1.0/(gamma*dt) : 1.0/((1.0 - gamma)*dt);

    /*--- We loop over faces in the range ---*/
    for(unsigned int face = face_range.first; face < face_range.second; ++face) {
      phi_proj_p.reinit(face);
      phi_proj_p.gather_evaluate(src[0], EvaluationFlags::values);
      phi_proj_m.reinit(face);
      phi_proj_m.gather_evaluate(src[0], EvaluationFlags::values);
      phi_p.reinit(face);
      phi_m.reinit(face);

      /*--- Now we loop over all the quadrature points to compute the integrals ---*/
      for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
        const auto& n_plus           = phi_p.get_normal_vector(q);
        const auto& avg_u_star_star  = 0.5*(phi_proj_p.get_value(q) + phi_proj_m.get_value(q));

        phi_p.submit_value(-coeff*scalar_product(avg_u_star_star, n_plus), q);
        phi_m.submit_value(coeff*scalar_product(avg_u_star_star, n_plus), q);
      }
      phi_p.integrate_scatter(EvaluationFlags::values, dst);
      phi_m.integrate_scatter(EvaluationFlags::values, dst);
    }
  }


  // The following function assembles rhs boundary term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_rhs_boundary_term_pressure(const MatrixFree<dim, Number>&               data,
                                      Vec&                                         dst,
                                      const std::vector<Vec>&                      src,
                                      const std::pair<unsigned int, unsigned int>& face_range) const {
    /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
    FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number>   phi(data, true, 1, 1);
    FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_p, dim, Number> phi_proj(data, true, 0, 1);

    const double coeff = (TR_BDF2_stage == 1) ? 1.0/(gamma*dt) : 1.0/((1.0 - gamma)*dt);

    /*--- We loop over faces in the range ---*/
    for(unsigned int face = face_range.first; face < face_range.second; ++face) {
      phi_proj.reinit(face);
      phi_proj.gather_evaluate(src[0], EvaluationFlags::values);
      phi.reinit(face);

      /*--- Now we loop over all the quadrature points to compute the integrals ---*/
      for(unsigned int q = 0; q < phi.n_q_points; ++q) {
        const auto& n_plus = phi.get_normal_vector(q);

        phi.submit_value(-coeff*scalar_product(phi_proj.get_value(q), n_plus), q);
      }
      phi.integrate_scatter(EvaluationFlags::values, dst);
    }
  }


  // Put together all the previous steps for pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  vmult_rhs_pressure(Vec& dst, const std::vector<Vec>& src) const {
    for(auto& vec : src)
      vec.update_ghost_values();

    this->data->loop(&NavierStokesProjectionOperator::assemble_rhs_cell_term_pressure,
                     &NavierStokesProjectionOperator::assemble_rhs_face_term_pressure,
                     &NavierStokesProjectionOperator::assemble_rhs_boundary_term_pressure,
                     this, dst, src, true,
                     MatrixFree<dim, Number>::DataAccessOnFaces::unspecified,
                     MatrixFree<dim, Number>::DataAccessOnFaces::unspecified);
  }


  // Now we need to build the 'matrices', i.e. the bilinear forms. We start by
  // assembling the cell term for the velocity
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_cell_term_velocity(const MatrixFree<dim, Number>&               data,
                              Vec&                                         dst,
                              const Vec&                                   src,
                              const std::pair<unsigned int, unsigned int>& cell_range) const {
    if(TR_BDF2_stage == 1) {
      /*--- We first start by declaring the suitable instances to read already available quantities. Moreover 'phi' in
      this case serves for a bilinear form and so it will not used only to submit but also to read the src ---*/
      FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0),
                                                                   phi_old_extr(data, 0);
      FEEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, 2),
                                                                   phi_y_plus(data, 2);

            /*--- We loop over all cells in the range ---*/
      for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
        phi.reinit(cell);
        phi.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_extr.reinit(cell);
        phi_old_extr.gather_evaluate(u_extr, EvaluationFlags::values);

        phi_deltas.reinit(cell);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(cell);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        /*--- Now we loop over all quadrature points ---*/
        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          const auto& u_int                = phi.get_value(q);
          const auto& grad_u_int           = 2.0 * phi.get_symmetric_gradient(q);
          const auto& u_n_gamma_ov_2       = phi_old_extr.get_value(q);
          const auto& tensor_product_u_int = outer_product(u_int, u_n_gamma_ov_2);

          const auto& dx                   = phi_deltas.get_value(q);

          phi.submit_value(1.0/(gamma*dt)*u_int, q);

          phi.submit_gradient(-a22*tensor_product_u_int + a22*viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re)*grad_u_int, q);
        }
        phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
    else {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0),
                                                                   phi_int_extr(data, 0);
      FEEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, 2),
                                                                   phi_y_plus(data, 2);

      /*--- We loop over all cells in the range ---*/
      for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
        phi.reinit(cell);
        phi.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_int_extr.reinit(cell);
        phi_int_extr.gather_evaluate(u_extr, EvaluationFlags::values);

        phi_deltas.reinit(cell);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(cell);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        /*--- Now we loop over all quadrature points ---*/
        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          const auto& u_curr                   = phi.get_value(q);
          const auto& grad_u_curr              = 2.0 * phi.get_symmetric_gradient(q);
          const auto& u_n1_int                 = phi_int_extr.get_value(q);
          const auto& tensor_product_u_curr    = outer_product(u_curr, u_n1_int);

          const auto& dx                       = phi_deltas.get_value(q);

          phi.submit_value(1.0/((1.0 - gamma)*dt)*u_curr, q);
          phi.submit_gradient(-a33*tensor_product_u_curr + a33*viscosity.value(phi_y_plus.get_value(q), grad_u_curr, dx, Re)*grad_u_curr, q);
        }
        phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
  }


  // The following function assembles face term for the velocity
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_face_term_velocity(const MatrixFree<dim, Number>&               data,
                              Vec&                                         dst,
                              const Vec&                                   src,
                              const std::pair<unsigned int, unsigned int>& face_range) const {
    if(TR_BDF2_stage == 1) {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi_p(data, true, 0),
                                                                       phi_m(data, false, 0),
                                                                       phi_old_extr_p(data, true, 0),
                                                                       phi_old_extr_m(data, false, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas_p(data, true, 2),
                                                                       phi_deltas_m(data, false, 2),
                                                                       phi_y_plus_p(data, true, 2),
                                                                       phi_y_plus_m(data, false, 2);

      /*--- We loop over all faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_p.reinit(face);
        phi_p.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_extr_p.reinit(face);
        phi_old_extr_p.gather_evaluate(u_extr, EvaluationFlags::values);
        phi_old_extr_m.reinit(face);
        phi_old_extr_m.gather_evaluate(u_extr, EvaluationFlags::values);

        phi_deltas_p.reinit(face);
        phi_deltas_p.gather_evaluate(deltas, EvaluationFlags::values);
        phi_deltas_m.reinit(face);
        phi_deltas_m.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus_p.reinit(face);
        phi_y_plus_p.gather_evaluate(y_plus, EvaluationFlags::values);
        phi_y_plus_m.reinit(face);
        phi_y_plus_m.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto coef_jump = C_u*0.5*(std::abs((phi_p.get_normal_vector(0)*phi_p.inverse_jacobian(0))[dim - 1]) +
                                        std::abs((phi_m.get_normal_vector(0)*phi_m.inverse_jacobian(0))[dim - 1]));

        /*--- Now we loop over all quadrature points ---*/
        for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
          const auto& n_plus                   = phi_p.get_normal_vector(q);

          const auto& dx_p                     = phi_deltas_p.get_value(q);
          const auto& dx_m                     = phi_deltas_m.get_value(q);

          const auto& avg_visc_grad_u_int      = viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                                 phi_p.get_symmetric_gradient(q) +
                                                 viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                                 phi_m.get_symmetric_gradient(q);
          const auto& avg_visc                 = 0.5 * (viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re) +
                                                        viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re));

          const auto& jump_u_int               = phi_p.get_value(q) - phi_m.get_value(q);
          const auto& avg_tensor_product_u_int = 0.5*(outer_product(phi_p.get_value(q), phi_old_extr_p.get_value(q)) +
                                                      outer_product(phi_m.get_value(q), phi_old_extr_m.get_value(q)));
          const auto  lambda                   = std::max(std::abs(scalar_product(phi_old_extr_p.get_value(q), n_plus)),
                                                          std::abs(scalar_product(phi_old_extr_m.get_value(q), n_plus)));

          phi_p.submit_value(a22*(-avg_visc_grad_u_int*n_plus + avg_visc * coef_jump*jump_u_int) +
                             a22*avg_tensor_product_u_int*n_plus + 0.5*a22*lambda*jump_u_int, q);
          phi_m.submit_value(-a22*(-avg_visc_grad_u_int*n_plus +  avg_visc * coef_jump*jump_u_int) -
                              a22*avg_tensor_product_u_int*n_plus - 0.5*a22*lambda*jump_u_int, q);

          phi_p.submit_gradient(-theta_v*a22*viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                 0.5*(outer_product(jump_u_int,n_plus) + outer_product(n_plus,jump_u_int)), q);
          phi_m.submit_gradient(-theta_v*a22*viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                 0.5*(outer_product(jump_u_int,n_plus) + outer_product(n_plus,jump_u_int)), q);

        }
        phi_p.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
        phi_m.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
    else {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi_p(data, true, 0),
                                                                       phi_m(data, false, 0),
                                                                       phi_extr_p(data, true, 0),
                                                                       phi_extr_m(data, false, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas_p(data, true, 2),
                                                                       phi_deltas_m(data, false, 2),
                                                                       phi_y_plus_p(data, true, 2),
                                                                       phi_y_plus_m(data, false, 2);

      /*--- We loop over all faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_p.reinit(face);
        phi_p.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_extr_p.reinit(face);
        phi_extr_p.gather_evaluate(u_extr, EvaluationFlags::values);
        phi_extr_m.reinit(face);
        phi_extr_m.gather_evaluate(u_extr, EvaluationFlags::values);

        phi_deltas_p.reinit(face);
        phi_deltas_p.gather_evaluate(deltas, EvaluationFlags::values);
        phi_deltas_m.reinit(face);
        phi_deltas_m.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus_p.reinit(face);
        phi_y_plus_p.gather_evaluate(y_plus, EvaluationFlags::values);
        phi_y_plus_m.reinit(face);
        phi_y_plus_m.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto coef_jump = C_u*0.5*(std::abs((phi_p.get_normal_vector(0)*phi_p.inverse_jacobian(0))[dim - 1]) +
                                        std::abs((phi_m.get_normal_vector(0)*phi_m.inverse_jacobian(0))[dim - 1]));

        /*--- Now we loop over all quadrature points ---*/
        for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
          const auto& n_plus               = phi_p.get_normal_vector(q);

          const auto& dx_p                 = phi_deltas_p.get_value(q);
          const auto& dx_m                 = phi_deltas_m.get_value(q);

          const auto& avg_visc_grad_u      = viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                             phi_p.get_symmetric_gradient(q) +
                                             viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                             phi_m.get_symmetric_gradient(q);
          const auto& avg_visc             = 0.5 * (viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re) +
                                                    viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re));

          const auto& jump_u               = phi_p.get_value(q) - phi_m.get_value(q);
          const auto& avg_tensor_product_u = 0.5*(outer_product(phi_p.get_value(q), phi_extr_p.get_value(q)) +
                                                  outer_product(phi_m.get_value(q), phi_extr_m.get_value(q)));
          const auto& lambda               = std::max(std::abs(scalar_product(phi_extr_p.get_value(q), n_plus)),
                                                      std::abs(scalar_product(phi_extr_m.get_value(q), n_plus)));

          phi_p.submit_value(a33*(-avg_visc_grad_u*n_plus + avg_visc*coef_jump*jump_u) +
                             a33*avg_tensor_product_u*n_plus + 0.5*a33*lambda*jump_u, q);
          phi_m.submit_value(-a33*(-avg_visc_grad_u*n_plus + avg_visc*coef_jump*jump_u) -
                              a33*avg_tensor_product_u*n_plus - 0.5*a33*lambda*jump_u, q);

          phi_p.submit_gradient(-theta_v*a33*viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                 0.5*(outer_product(jump_u,n_plus) + outer_product(n_plus,jump_u)), q);
          phi_m.submit_gradient(-theta_v*a33*viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                 0.5*(outer_product(jump_u,n_plus) + outer_product(n_plus,jump_u)), q);

        }
        phi_p.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
        phi_m.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
  }


  // The following function assembles boundary term for the velocity
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_boundary_term_velocity(const MatrixFree<dim, Number>&               data,
                                  Vec&                                         dst,
                                  const Vec&                                   src,
                                  const std::pair<unsigned int, unsigned int>& face_range) const {
    if(TR_BDF2_stage == 1) {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, true, 0),
                                                                       phi_old_extr(data, true, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, true, 2),
                                                                       phi_y_plus(data, true, 2);

      /*--- We loop over all faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi.reinit(face);
        phi.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_old_extr.reinit(face);
        phi_old_extr.gather_evaluate(u_extr, EvaluationFlags::values);

        phi_deltas.reinit(face);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(face);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto boundary_id = data.get_boundary_id(face);
        const auto coef_jump   = C_u*std::abs((phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]);

        /*--- The application of the mirror principle is not so trivial because we have a Dirichlet condition
              on a single component for the outflow; so we distinguish the two cases ---*/
        if(boundary_id != 1 && (boundary_id != 3 || no_slip)) {
          const double coef_trasp = 0.0;

          /*--- Now we loop over all quadrature points ---*/
          for(unsigned int q = 0; q < phi.n_q_points; ++q) {
            const auto& n_plus               = phi.get_normal_vector(q);
            const auto& grad_u_int           = 2.0 * phi.get_symmetric_gradient(q);
            const auto& u_int                = phi.get_value(q);
            const auto& tensor_product_u_int = outer_product(phi.get_value(q), phi_old_extr.get_value(q));
            const auto& lambda               = std::abs(scalar_product(phi_old_extr.get_value(q), n_plus));

            const auto& dx                   = phi_deltas.get_value(q);

            phi.submit_value(a22*viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re) * (-grad_u_int*n_plus + 2.0*coef_jump*u_int) +
                             a22*coef_trasp*tensor_product_u_int*n_plus + a22*lambda*u_int, q);

            phi.submit_gradient(-theta_v*a22*viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re)*
                                 (outer_product(u_int, n_plus) + outer_product(n_plus, u_int)), q);
          }
          phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
        }
        else {
          /*--- Now we loop over all quadrature points ---*/
          for(unsigned int q = 0; q < phi.n_q_points; ++q) {
            const auto& n_plus               = phi.get_normal_vector(q);
            const auto& grad_u_int           = 2.0 * phi.get_symmetric_gradient(q);
            const auto& u_int                = phi.get_value(q);
            const auto& lambda               = std::abs(scalar_product(phi_old_extr.get_value(q), n_plus));
            const auto& dx = phi_deltas.get_value(q);

            const auto& point_vectorized     = phi.quadrature_point(q);
            auto u_int_m                     = u_int;
            auto grad_u_int_m                = grad_u_int;
            for(unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
              Point<dim> point;
              for(unsigned int d = 0; d < dim; ++d)
                point[d] = point_vectorized[d][v];

              u_int_m[1][v] = -u_int_m[1][v];

              grad_u_int_m[0][0][v] = -grad_u_int_m[0][0][v];
              grad_u_int_m[0][1][v] = -grad_u_int_m[0][1][v];
            }

            const auto& visc = viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re);

            phi.submit_value(a22*(-(0.5*(grad_u_int*visc + grad_u_int_m*viscosity.value(phi_y_plus.get_value(q), grad_u_int_m, dx, Re)))*n_plus +
                             visc*coef_jump*(u_int - u_int_m)) + a22*outer_product(0.5*(u_int + u_int_m), phi_old_extr.get_value(q))*n_plus +
                             a22*0.5*lambda*(u_int - u_int_m), q);
            phi.submit_gradient(-theta_v*a22*visc*(outer_product(u_int - u_int_m, n_plus) + outer_product(n_plus, u_int - u_int_m)), q);
          }
          phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
        }
      }
    }
    else {
      /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, true, 0),
                                                                       phi_extr(data, true, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, true, 2),
                                                                       phi_y_plus(data, true, 2);


      /*--- We loop over all faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi.reinit(face);
        phi.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
        phi_extr.reinit(face);
        phi_extr.gather_evaluate(u_extr, EvaluationFlags::values);

        phi_deltas.reinit(face);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(face);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto boundary_id = data.get_boundary_id(face);
        const auto coef_jump   = C_u*std::abs((phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]);

        if(boundary_id != 1 && (boundary_id != 3 || no_slip)) {
          const double coef_trasp = 0.0;

          /*--- Now we loop over all quadrature points ---*/
          for(unsigned int q = 0; q < phi.n_q_points; ++q) {
            const auto& n_plus           = phi.get_normal_vector(q);
            const auto& grad_u           = 2. * phi.get_symmetric_gradient(q);
            const auto& u                = phi.get_value(q);
            const auto& tensor_product_u = outer_product(phi.get_value(q), phi_extr.get_value(q));
            const auto& lambda           = std::abs(scalar_product(phi_extr.get_value(q), n_plus));

            const auto& dx               = phi_deltas.get_value(q);

            phi.submit_value(a33*viscosity.value(phi_y_plus.get_value(q), grad_u, dx, Re)*(-grad_u*n_plus + 2.0*coef_jump*u) +
                             a33*coef_trasp*tensor_product_u*n_plus + a33*lambda*u, q);
            phi.submit_gradient(-theta_v*a33*viscosity.value(phi_y_plus.get_value(q), grad_u, dx, Re)*
                                 (outer_product(u, n_plus) + outer_product(n_plus, u)), q);

          }
          phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
        }
        else {
          /*--- Now we loop over all quadrature points ---*/
          for(unsigned int q = 0; q < phi.n_q_points; ++q) {
            const auto& n_plus           = phi.get_normal_vector(q);
            const auto& grad_u           = 2.0 * phi.get_symmetric_gradient(q);
            const auto& u                = phi.get_value(q);
            const auto& lambda           = std::abs(scalar_product(phi_extr.get_value(q), n_plus));

            const auto& dx               = phi_deltas.get_value(q);
            const auto& point_vectorized = phi.quadrature_point(q);
            auto u_m                     = u;
            auto grad_u_m                = grad_u;
            for(unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
              Point<dim> point;
              for(unsigned int d = 0; d < dim; ++d)
                point[d] = point_vectorized[d][v];

              u_m[1][v] = -u_m[1][v];

              grad_u_m[0][0][v] = -grad_u_m[0][0][v];
              grad_u_m[0][1][v] = -grad_u_m[0][1][v];
            }
            const auto& visc = viscosity.value(phi_y_plus.get_value(q), grad_u, dx, Re);

            phi.submit_value(a33*(-(0.5*(visc*grad_u + viscosity.value(phi_y_plus.get_value(q), grad_u_m, dx, Re)*grad_u_m))*n_plus +
                             visc*coef_jump*(u - u_m)) + a33*outer_product(0.5*(u + u_m), phi_extr.get_value(q))*n_plus +
                             a33*0.5*lambda*(u - u_m), q);
            phi.submit_gradient(-theta_v*a33*visc*(outer_product(u - u_m, n_plus) + outer_product(n_plus, u - u_m)), q);
          }
          phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
        }
      }
    }
  }


  // Next, we focus on 'matrices' to compute the pressure. We first assemble cell term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_cell_term_pressure(const MatrixFree<dim, Number>&               data,
                              Vec&                                         dst,
                              const Vec&                                   src,
                              const std::pair<unsigned int, unsigned int>& cell_range) const {
    /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
    FEEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number> phi(data, 1, 1);

    const double coeff = (TR_BDF2_stage == 1) ? 1.0e6*gamma*dt*gamma*dt : 1.0e6*(1.0 - gamma)*dt*(1.0 - gamma)*dt;

    /*--- Loop over all cells in the range ---*/
    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
      phi.reinit(cell);
      phi.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);

      /*--- Now we loop over all quadrature points ---*/
      for(unsigned int q = 0; q < phi.n_q_points; ++q) {
        phi.submit_gradient(phi.get_gradient(q), q);
        phi.submit_value(1.0/coeff*phi.get_value(q), q);
      }

      phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
    }
  }


  // The following function assembles face term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_face_term_pressure(const MatrixFree<dim, Number>&               data,
                              Vec&                                         dst,
                              const Vec&                                   src,
                              const std::pair<unsigned int, unsigned int>& face_range) const {
    /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
    FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number> phi_p(data, true, 1, 1),
                                                                   phi_m(data, false, 1, 1);

    /*--- Loop over all faces in the range ---*/
    for(unsigned int face = face_range.first; face < face_range.second; ++face) {
      phi_p.reinit(face);
      phi_p.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);
      phi_m.reinit(face);
      phi_m.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);

      const auto coef_jump = C_p*0.5*(std::abs((phi_p.get_normal_vector(0)*phi_p.inverse_jacobian(0))[dim - 1]) +
                                      std::abs((phi_m.get_normal_vector(0)*phi_m.inverse_jacobian(0))[dim - 1]));

      /*--- Loop over quadrature points ---*/
      for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
        const auto& n_plus        = phi_p.get_normal_vector(q);

        const auto& avg_grad_pres = 0.5*(phi_p.get_gradient(q) + phi_m.get_gradient(q));
        const auto& jump_pres     = phi_p.get_value(q) - phi_m.get_value(q);

        phi_p.submit_value(-scalar_product(avg_grad_pres, n_plus) + coef_jump*jump_pres, q);
        phi_m.submit_value(scalar_product(avg_grad_pres, n_plus) - coef_jump*jump_pres, q);
        phi_p.submit_gradient(-theta_p*0.5*jump_pres*n_plus, q);
        phi_m.submit_gradient(-theta_p*0.5*jump_pres*n_plus, q);
      }
      phi_p.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      phi_m.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
    }
  }


  // The following function assembles boundary term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_boundary_term_pressure(const MatrixFree<dim, Number>&               data,
                                  Vec&                                         dst,
                                  const Vec&                                   src,
                                  const std::pair<unsigned int, unsigned int>& face_range) const {
    FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number> phi(data, true, 1, 1);

    for(unsigned int face = face_range.first; face < face_range.second; ++face) {
      const auto boundary_id = data.get_boundary_id(face);

      if(boundary_id == 1|| (!no_slip && boundary_id == 3)) {
        phi.reinit(face);
        phi.gather_evaluate(src, true, true);

        const auto coef_jump = C_p*std::abs((phi.get_normal_vector(0)*phi.inverse_jacobian(0))[dim - 1]);

        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          const auto& n_plus    = phi.get_normal_vector(q);

          const auto& grad_pres = phi.get_gradient(q);
          const auto& pres      = phi.get_value(q);

          phi.submit_value(-scalar_product(grad_pres, n_plus) + coef_jump*pres , q);
          phi.submit_normal_derivative(-theta_p*pres, q);
        }
        phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
      }
    }
  }


  // Before coding the 'apply_add' function, which is the one that will perform the loop, we focus on
  // the linear system that arises to project the gradient of the pressure into the velocity space.
  // The following function assembles rhs cell term for the projection of gradient of pressure. Since no
  // integration by parts is performed, only a cell term contribution is present.
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_rhs_cell_term_projection_grad_p(const MatrixFree<dim, Number>&               data,
                                           Vec&                                         dst,
                                           const Vec&                                   src,
                                           const std::pair<unsigned int, unsigned int>& cell_range) const {
    /*--- We first start by declaring the suitable instances to read already available quantities. ---*/
    FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0);
    FEEvaluation<dim, fe_degree_p, n_q_points_1d_v, 1, Number>   phi_pres(data, 1);

    /*--- Loop over all cells in the range ---*/
    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
      phi_pres.reinit(cell);
      phi_pres.gather_evaluate(src, EvaluationFlags::gradients);
      phi.reinit(cell);

      /*--- Loop over quadrature points ---*/
      for(unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_value(phi_pres.get_gradient(q), q);

      phi.integrate_scatter(EvaluationFlags::values, dst);
    }
  }


  // Put together all the previous steps for porjection of pressure gradient. Here we loop only over cells
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  vmult_grad_p_projection(Vec& dst, const Vec& src) const {
    this->data->cell_loop(&NavierStokesProjectionOperator::assemble_rhs_cell_term_projection_grad_p,
                          this, dst, src, true);
  }


  // Assemble now cell term for the projection of gradient of pressure. This is nothing but a mass matrix
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_cell_term_projection_grad_p(const MatrixFree<dim, Number>&               data,
                                       Vec&                                         dst,
                                       const Vec&                                   src,
                                       const std::pair<unsigned int, unsigned int>& cell_range) const {
    FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0);

    /*--- Loop over all cells in the range ---*/
    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
      phi.reinit(cell);
      phi.gather_evaluate(src, EvaluationFlags::values);

      /*--- Loop over quadrature points ---*/
      for(unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_value(phi.get_value(q), q);

      phi.integrate_scatter(EvaluationFlags::values, dst);
    }
  }


  // Put together all previous steps. This is the overriden function that effectively performs the
  // matrix-vector multiplication.
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  apply_add(Vec& dst, const Vec& src) const {
    if(NS_stage == 1) {
      this->data->loop(&NavierStokesProjectionOperator::assemble_cell_term_velocity,
                       &NavierStokesProjectionOperator::assemble_face_term_velocity,
                       &NavierStokesProjectionOperator::assemble_boundary_term_velocity,
                       this, dst, src, false,
                       MatrixFree<dim, Number>::DataAccessOnFaces::unspecified,
                       MatrixFree<dim, Number>::DataAccessOnFaces::unspecified);
    }
    else if(NS_stage == 2) {
      this->data->loop(&NavierStokesProjectionOperator::assemble_cell_term_pressure,
                       &NavierStokesProjectionOperator::assemble_face_term_pressure,
                       &NavierStokesProjectionOperator::assemble_boundary_term_pressure,
                       this, dst, src, false,
                       MatrixFree<dim, Number>::DataAccessOnFaces::unspecified,
                       MatrixFree<dim, Number>::DataAccessOnFaces::unspecified);
    }
    else if(NS_stage == 3) {
      this->data->cell_loop(&NavierStokesProjectionOperator::assemble_cell_term_projection_grad_p,
                            this, dst, src, false); /*--- Since we have only a cell term contribution, we use cell_loop ---*/
    }
    else
      Assert(false, ExcNotImplemented());
  }


  // Finally, we focus on computing the diagonal for preconditioners and we start by assembling
  // the diagonal cell term for the velocity. Since we do not have access to the entries of the matrix,
  // in order to compute the element i, we test the matrix against a vector which is equal to 1 in position i and 0 elsewhere.
  // This is why 'src' will result as unused.
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_diagonal_cell_term_velocity(const MatrixFree<dim, Number>&               data,
                                       Vec&                                         dst,
                                       const unsigned int&                          ,
                                       const std::pair<unsigned int, unsigned int>& cell_range) const {
    if(TR_BDF2_stage == 1) {
      FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0),
                                                                   phi_old_extr(data, 0);
      FEEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, 2),
                                                                   phi_y_plus(data, 2);

      AlignedVector<Tensor<1, dim, VectorizedArray<Number>>> diagonal(phi.dofs_per_component);
      /*--- Build a vector of ones to be tested (here we will see the velocity as a whole vector, since
                                                 dof_handler_velocity is vectorial and so the dof values are vectors). ---*/
      Tensor<1, dim, VectorizedArray<Number>> tmp;
      for(unsigned int d = 0; d < dim; ++d)
        tmp[d] = make_vectorized_array<Number>(1.0);

      /*--- Loop over cells in the range ---*/
      for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
        phi_old_extr.reinit(cell);
        phi_old_extr.gather_evaluate(u_extr, true, false);
        phi.reinit(cell);

        phi_deltas.reinit(cell);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(cell);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        /*--- Loop over dofs ---*/
        for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
          for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
            phi.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j); /*--- Set all dofs to zero ---*/
          phi.submit_dof_value(tmp, i); /*--- Set dof i equal to one ---*/
          phi.evaluate(true, true);

          /*--- Loop over quadrature points ---*/
          for(unsigned int q = 0; q < phi.n_q_points; ++q) {
            const auto& u_int                = phi.get_value(q);
            const auto& grad_u_int           = 2.0 * phi.get_symmetric_gradient(q);
            const auto& u_n_gamma_ov_2       = phi_old_extr.get_value(q);
            const auto& tensor_product_u_int = outer_product(u_int, u_n_gamma_ov_2);

            const auto& dx                   = phi_deltas.get_value(q);

            phi.submit_value(1.0/(gamma*dt)*u_int, q);
            phi.submit_gradient(-a22*tensor_product_u_int + a22*viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re)*grad_u_int, q);
          }
          phi.integrate(true, true);
          diagonal[i] = phi.get_dof_value(i);
        }
        for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
          phi.submit_dof_value(diagonal[i], i);
        phi.distribute_local_to_global(dst);
      }
    }
    else {
      FEEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, 0),
                                                                   phi_int_extr(data, 0);
      FEEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, 2),
                                                                   phi_y_plus(data, 2);

      AlignedVector<Tensor<1, dim, VectorizedArray<Number>>> diagonal(phi.dofs_per_component);
      Tensor<1, dim, VectorizedArray<Number>> tmp;
      for(unsigned int d = 0; d < dim; ++d)
        tmp[d] = make_vectorized_array<Number>(1.0);

      /*--- Loop over cells in the range ---*/
      for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
        phi_int_extr.reinit(cell);
        phi_int_extr.gather_evaluate(u_extr, true, false);
        phi.reinit(cell);

        phi_deltas.reinit(cell);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(cell);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        /*--- Loop over dofs ---*/
        for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
          for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
            phi.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
          phi.submit_dof_value(tmp, i);
          phi.evaluate(true, true);

          /*--- Loop over quadrature points ---*/
          for(unsigned int q = 0; q < phi.n_q_points; ++q) {
            const auto& u_curr                   = phi.get_value(q);
            const auto& grad_u_curr              = 2.0 * phi.get_symmetric_gradient(q);
            const auto& u_n1_int                 = phi_int_extr.get_value(q);
            const auto& tensor_product_u_curr    = outer_product(u_curr, u_n1_int);

            const auto& dx                       = phi_deltas.get_value(q);

            phi.submit_value(1.0/((1.0 - gamma)*dt)*u_curr, q);
            phi.submit_gradient(-a33*tensor_product_u_curr + a33*viscosity.value(phi_y_plus.get_value(q), grad_u_curr, dx, Re)*grad_u_curr, q);
          }
          phi.integrate(true, true);
          diagonal[i] = phi.get_dof_value(i);
        }
        for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
          phi.submit_dof_value(diagonal[i], i);
        phi.distribute_local_to_global(dst);
      }
    }
  }


  // The following function assembles diagonal face term for the velocity
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_diagonal_face_term_velocity(const MatrixFree<dim, Number>&               data,
                                       Vec&                                         dst,
                                       const unsigned int&                          ,
                                       const std::pair<unsigned int, unsigned int>& face_range) const {
    if(TR_BDF2_stage == 1) {
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi_p(data, true, 0),
                                                                       phi_m(data, false, 0),
                                                                       phi_old_extr_p(data, true, 0),
                                                                       phi_old_extr_m(data, false, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas_p(data, true, 2),
                                                                       phi_deltas_m(data, false, 2),
                                                                       phi_y_plus_p(data, true, 2),
                                                                       phi_y_plus_m(data, false, 2);

      AssertDimension(phi_p.dofs_per_component, phi_m.dofs_per_component); /*--- We just assert for safety that dimension match,
                                                                                in the sense that we have selected the proper
                                                                                space ---*/
      AlignedVector<Tensor<1, dim, VectorizedArray<Number>>> diagonal_p(phi_p.dofs_per_component),
                                                             diagonal_m(phi_m.dofs_per_component);

      Tensor<1, dim, VectorizedArray<Number>> tmp;
      for(unsigned int d = 0; d < dim; ++d)
        tmp[d] = make_vectorized_array<Number>(1.0); /*--- We build the usal vector of ones that we will use as dof value ---*/

      /*--- Now we loop over faces ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_old_extr_p.reinit(face);
        phi_old_extr_p.gather_evaluate(u_extr, true, false);
        phi_old_extr_m.reinit(face);
        phi_old_extr_m.gather_evaluate(u_extr, true, false);
        phi_p.reinit(face);
        phi_m.reinit(face);

        phi_deltas_p.reinit(face);
        phi_deltas_p.gather_evaluate(deltas, EvaluationFlags::values);
        phi_deltas_m.reinit(face);
        phi_deltas_m.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus_p.reinit(face);
        phi_y_plus_p.gather_evaluate(y_plus, EvaluationFlags::values);
        phi_y_plus_m.reinit(face);
        phi_y_plus_m.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto coef_jump = C_u*0.5*(std::abs((phi_p.get_normal_vector(0)*phi_p.inverse_jacobian(0))[dim - 1]) +
                                        std::abs((phi_m.get_normal_vector(0)*phi_m.inverse_jacobian(0))[dim - 1]));

        /*--- Loop over dofs. We will set all equal to zero apart from the current one ---*/
        for(unsigned int i = 0; i < phi_p.dofs_per_component; ++i) {
          for(unsigned int j = 0; j < phi_p.dofs_per_component; ++j) {
            phi_p.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
            phi_m.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
          }
          phi_p.submit_dof_value(tmp, i);
          phi_p.evaluate(true, true);
          phi_m.submit_dof_value(tmp, i);
          phi_m.evaluate(true, true);

          /*--- Loop over quadrature points to compute the integral ---*/
          for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
            const auto& n_plus                   = phi_p.get_normal_vector(q);

            const auto& dx_p                     = phi_deltas_p.get_value(q);
            const auto& dx_m                     = phi_deltas_m.get_value(q);

            const auto& avg_visc_grad_u_int      = viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                                   phi_p.get_symmetric_gradient(q) +
                                                   viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                                   phi_m.get_symmetric_gradient(q);
            const auto& avg_visc                 = 0.5 * (viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re) +
                                                          viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re));

            const auto& jump_u_int               = phi_p.get_value(q) - phi_m.get_value(q);
            const auto& avg_tensor_product_u_int = 0.5*(outer_product(phi_p.get_value(q), phi_old_extr_p.get_value(q)) +
                                                        outer_product(phi_m.get_value(q), phi_old_extr_m.get_value(q)));
            const auto& lambda                   = std::max(std::abs(scalar_product(phi_old_extr_p.get_value(q), n_plus)),
                                                            std::abs(scalar_product(phi_old_extr_m.get_value(q), n_plus)));


          phi_p.submit_value(a22*(-avg_visc_grad_u_int*n_plus + avg_visc * coef_jump*jump_u_int) +
                             a22*avg_tensor_product_u_int*n_plus + 0.5*a22*lambda*jump_u_int, q);
          phi_m.submit_value(-a22*(-avg_visc_grad_u_int*n_plus +  avg_visc * coef_jump*jump_u_int) -
                              a22*avg_tensor_product_u_int*n_plus - 0.5*a22*lambda*jump_u_int, q);

          phi_p.submit_gradient(-theta_v*a22*viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                 0.5*(outer_product(jump_u_int,n_plus) + outer_product(n_plus,jump_u_int)), q);
          phi_m.submit_gradient(-theta_v*a22*viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                 0.5*(outer_product(jump_u_int,n_plus) + outer_product(n_plus,jump_u_int)), q);
          }
          phi_p.integrate(true, true);
          diagonal_p[i] = phi_p.get_dof_value(i);
          phi_m.integrate(true, true);
          diagonal_m[i] = phi_m.get_dof_value(i);
        }
        for(unsigned int i = 0; i < phi_p.dofs_per_component; ++i) {
          phi_p.submit_dof_value(diagonal_p[i], i);
          phi_m.submit_dof_value(diagonal_m[i], i);
        }
        phi_p.distribute_local_to_global(dst);
        phi_m.distribute_local_to_global(dst);
      }
    }
    else {
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi_p(data, true, 0),
                                                                       phi_m(data, false, 0),
                                                                       phi_extr_p(data, true, 0),
                                                                       phi_extr_m(data, false, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas_p(data, true, 2),
                                                                       phi_deltas_m(data, false, 2),
                                                                       phi_y_plus_p(data, true, 2),
                                                                       phi_y_plus_m(data, false, 2);

      AssertDimension(phi_p.dofs_per_component, phi_m.dofs_per_component);
      AlignedVector<Tensor<1, dim, VectorizedArray<Number>>> diagonal_p(phi_p.dofs_per_component),
                                                             diagonal_m(phi_m.dofs_per_component);
      Tensor<1, dim, VectorizedArray<Number>> tmp;
      for(unsigned int d = 0; d < dim; ++d)
        tmp[d] = make_vectorized_array<Number>(1.0);

      /*--- Now we loop over faces ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_extr_p.reinit(face);
        phi_extr_p.gather_evaluate(u_extr, true, false);
        phi_extr_m.reinit(face);
        phi_extr_m.gather_evaluate(u_extr, true, false);
        phi_p.reinit(face);
        phi_m.reinit(face);

        phi_deltas_p.reinit(face);
        phi_deltas_p.gather_evaluate(deltas, EvaluationFlags::values);
        phi_deltas_m.reinit(face);
        phi_deltas_m.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus_p.reinit(face);
        phi_y_plus_p.gather_evaluate(y_plus, EvaluationFlags::values);
        phi_y_plus_m.reinit(face);
        phi_y_plus_m.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto coef_jump = C_u*0.5*(std::abs((phi_p.get_normal_vector(0)*phi_p.inverse_jacobian(0))[dim - 1]) +
                                        std::abs((phi_m.get_normal_vector(0)*phi_m.inverse_jacobian(0))[dim - 1]));

        /*--- Loop over dofs. We will set all equal to zero apart from the current one ---*/
        for(unsigned int i = 0; i < phi_p.dofs_per_component; ++i) {
          for(unsigned int j = 0; j < phi_p.dofs_per_component; ++j) {
            phi_p.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
            phi_m.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
          }
          phi_p.submit_dof_value(tmp, i);
          phi_p.evaluate(true, true);
          phi_m.submit_dof_value(tmp, i);
          phi_m.evaluate(true, true);

          /*--- Loop over quadrature points to compute the integral ---*/
          for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
            const auto& n_plus               = phi_p.get_normal_vector(q);

            const auto& dx_p                 = phi_deltas_p.get_value(q);
            const auto& dx_m                 = phi_deltas_m.get_value(q);

            const auto& avg_visc_grad_u      = viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                               phi_p.get_symmetric_gradient(q) +
                                               viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                               phi_m.get_symmetric_gradient(q);
            const auto& avg_visc = 0.5 * (viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re) +
                                          viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re));

            const auto& jump_u               = phi_p.get_value(q) - phi_m.get_value(q);
            const auto& avg_tensor_product_u = 0.5*(outer_product(phi_p.get_value(q), phi_extr_p.get_value(q)) +
                                                    outer_product(phi_m.get_value(q), phi_extr_m.get_value(q)));
            const auto& lambda               = std::max(std::abs(scalar_product(phi_extr_p.get_value(q), n_plus)),
                                                        std::abs(scalar_product(phi_extr_m.get_value(q), n_plus)));

            phi_p.submit_value(a33*(-avg_visc_grad_u*n_plus + avg_visc*coef_jump*jump_u) +
                               a33*avg_tensor_product_u*n_plus + 0.5*a33*lambda*jump_u, q);
            phi_m.submit_value(-a33*(-avg_visc_grad_u*n_plus + avg_visc*coef_jump*jump_u) -
                                a33*avg_tensor_product_u*n_plus - 0.5*a33*lambda*jump_u, q);

            phi_p.submit_gradient(-theta_v*a33*viscosity.value(phi_y_plus_p.get_value(q), phi_p.get_symmetric_gradient(q), dx_p, Re)*
                                   0.5*(outer_product(jump_u,n_plus) + outer_product(n_plus,jump_u)), q);
            phi_m.submit_gradient(-theta_v*a33*viscosity.value(phi_y_plus_m.get_value(q), phi_m.get_symmetric_gradient(q), dx_m, Re)*
                                   0.5*(outer_product(jump_u,n_plus) + outer_product(n_plus,jump_u)), q);
          }
          phi_p.integrate(true, true);
          diagonal_p[i] = phi_p.get_dof_value(i);
          phi_m.integrate(true, true);
          diagonal_m[i] = phi_m.get_dof_value(i);
        }
        for(unsigned int i = 0; i < phi_p.dofs_per_component; ++i) {
          phi_p.submit_dof_value(diagonal_p[i], i);
          phi_m.submit_dof_value(diagonal_m[i], i);
        }
        phi_p.distribute_local_to_global(dst);
        phi_m.distribute_local_to_global(dst);
      }
    }
  }


  // The following function assembles boundary term for the velocity
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_diagonal_boundary_term_velocity(const MatrixFree<dim, Number>&               data,
                                           Vec&                                         dst,
                                           const unsigned int&                          ,
                                           const std::pair<unsigned int, unsigned int>& face_range) const {
    if(TR_BDF2_stage == 1) {
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, true, 0),
                                                                       phi_old_extr(data, true, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, true, 2),
                                                                       phi_y_plus(data, true, 2);

      AlignedVector<Tensor<1, dim, VectorizedArray<Number>>> diagonal(phi.dofs_per_component);
      Tensor<1, dim, VectorizedArray<Number>> tmp;
      for(unsigned int d = 0; d < dim; ++d)
        tmp[d] = make_vectorized_array<Number>(1.0);

      /*--- Loop over all faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_old_extr.reinit(face);
        phi_old_extr.gather_evaluate(u_extr, true, false);
        phi.reinit(face);

        phi_deltas.reinit(face);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(face);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto boundary_id = data.get_boundary_id(face);
        const auto coef_jump   = C_u*std::abs((phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]);

        if(boundary_id != 1) {
          const double coef_trasp = 0.0;

          /*--- Loop over all dofs ---*/
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
            for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
              phi.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
            phi.submit_dof_value(tmp, i);
            phi.evaluate(true, true);

            /*--- Loop over quadrature points to compute the integral ---*/
            for(unsigned int q = 0; q < phi.n_q_points; ++q) {
              const auto& n_plus               = phi.get_normal_vector(q);
              const auto& grad_u_int           = 2.0 * phi.get_symmetric_gradient(q);
              const auto& u_int                = phi.get_value(q);
              const auto& tensor_product_u_int = outer_product(phi.get_value(q), phi_old_extr.get_value(q));
              const auto& lambda               = std::abs(scalar_product(phi_old_extr.get_value(q), n_plus));

              const auto& dx                   = phi_deltas.get_value(q);

              phi.submit_value(a22*viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re) * (-grad_u_int*n_plus + 2.0*coef_jump*u_int) +
                               a22*coef_trasp*tensor_product_u_int*n_plus + a22*lambda*u_int, q);

              phi.submit_gradient(-theta_v*a22*viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re)*
                                   (outer_product(u_int, n_plus)+ outer_product(n_plus, u_int)), q);
            }
            phi.integrate(true, true);
            diagonal[i] = phi.get_dof_value(i);
          }
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
            phi.submit_dof_value(diagonal[i], i);
          phi.distribute_local_to_global(dst);
        }
        else {
          /*--- Loop over all dofs ---*/
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
            for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
              phi.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
            phi.submit_dof_value(tmp, i);
            phi.evaluate(true, true);

            /*--- Loop over quadrature points to compute the integral ---*/
            for(unsigned int q = 0; q < phi.n_q_points; ++q) {
              const auto& n_plus               = phi.get_normal_vector(q);
              const auto& grad_u_int           = 2.0 * phi.get_symmetric_gradient(q);
              const auto& u_int                = phi.get_value(q);
              const auto& lambda               = std::abs(scalar_product(phi_old_extr.get_value(q), n_plus));

              const auto& dx                   = phi_deltas.get_value(q);

              const auto& point_vectorized     = phi.quadrature_point(q);
              auto u_int_m                     = u_int;
              auto grad_u_int_m                = grad_u_int;
              for(unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
                Point<dim> point;
                for(unsigned int d = 0; d < dim; ++d)
                  point[d] = point_vectorized[d][v];

                u_int_m[1][v] = -u_int_m[1][v];

                grad_u_int_m[0][0][v] = -grad_u_int_m[0][0][v];
                grad_u_int_m[0][1][v] = -grad_u_int_m[0][1][v];
              }

              const auto& visc = viscosity.value(phi_y_plus.get_value(q), grad_u_int, dx, Re);

              phi.submit_value(a22*(-(0.5*(grad_u_int*visc + grad_u_int_m*viscosity.value(phi_y_plus.get_value(q), grad_u_int_m, dx, Re)))*n_plus +
                              visc*coef_jump*(u_int - u_int_m)) + a22*outer_product(0.5*(u_int + u_int_m), phi_old_extr.get_value(q))*n_plus +
                              a22*0.5*lambda*(u_int - u_int_m), q);
              phi.submit_gradient(-theta_v*a22*visc*
                                   (outer_product(u_int - u_int_m, n_plus)+outer_product(n_plus, u_int - u_int_m)), q);
            }
            phi.integrate(true, true);
            diagonal[i] = phi.get_dof_value(i);
          }
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
            phi.submit_dof_value(diagonal[i], i);
          phi.distribute_local_to_global(dst);
        }
      }
    }
    else {
      FEFaceEvaluation<dim, fe_degree_v, n_q_points_1d_v, dim, Number> phi(data, true, 0),
                                                                       phi_extr(data, true, 0);
      FEFaceEvaluation<dim, 0, n_q_points_1d_v, 1, Number>             phi_deltas(data, true, 2),
                                                                       phi_y_plus(data, true, 2);

      AlignedVector<Tensor<1, dim, VectorizedArray<Number>>> diagonal(phi.dofs_per_component);
      Tensor<1, dim, VectorizedArray<Number>> tmp;
      for(unsigned int d = 0; d < dim; ++d)
        tmp[d] = make_vectorized_array<Number>(1.0);

      /*--- Loop over all faces in the range ---*/
      for(unsigned int face = face_range.first; face < face_range.second; ++face) {
        phi_extr.reinit(face);
        phi_extr.gather_evaluate(u_extr, true, false);
        phi.reinit(face);

        phi_deltas.reinit(face);
        phi_deltas.gather_evaluate(deltas, EvaluationFlags::values);

        phi_y_plus.reinit(face);
        phi_y_plus.gather_evaluate(y_plus, EvaluationFlags::values);

        const auto boundary_id = data.get_boundary_id(face);
        const auto coef_jump   = C_u*std::abs((phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]);

        if(boundary_id != 1) {
          const double coef_trasp = 0.0;

          /*--- Loop over all dofs ---*/
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
            for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
              phi.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
            phi.submit_dof_value(tmp, i);
            phi.evaluate(true, true);

            /*--- Loop over quadrature points to compute the integral ---*/
            for(unsigned int q = 0; q < phi.n_q_points; ++q) {
              const auto& n_plus           = phi.get_normal_vector(q);
              const auto& grad_u           = 2. * phi.get_symmetric_gradient(q);
              const auto& u                = phi.get_value(q);
              const auto& tensor_product_u = outer_product(phi.get_value(q), phi_extr.get_value(q));
              const auto& lambda           = std::abs(scalar_product(phi_extr.get_value(q), n_plus));

              const auto& dx               = phi_deltas.get_value(q);

              phi.submit_value(a33*viscosity.value(phi_y_plus.get_value(q), grad_u, dx, Re)*(-grad_u*n_plus + 2.0*coef_jump*u) +
                              a33*coef_trasp*tensor_product_u*n_plus + a33*lambda*u, q);
              phi.submit_gradient(-theta_v*a33*viscosity.value(phi_y_plus.get_value(q), grad_u, dx, Re)*
                                   (outer_product(u, n_plus) + outer_product(n_plus, u)), q);
            }
            phi.integrate(true, true);
            diagonal[i] = phi.get_dof_value(i);
          }
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
            phi.submit_dof_value(diagonal[i], i);
          phi.distribute_local_to_global(dst);
        }
        else {
          /*--- Loop over all dofs ---*/
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
            for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
              phi.submit_dof_value(Tensor<1, dim, VectorizedArray<Number>>(), j);
            phi.submit_dof_value(tmp, i);
            phi.evaluate(true, true);

            /*--- Loop over quadrature points to compute the integral ---*/
            for(unsigned int q = 0; q < phi.n_q_points; ++q) {
              const auto& n_plus           = phi.get_normal_vector(q);
              const auto& grad_u           = 2.0 * phi.get_symmetric_gradient(q);
              const auto& u                = phi.get_value(q);
              const auto& lambda           = std::abs(scalar_product(phi_extr.get_value(q), n_plus));

              const auto& dx = phi_deltas.get_value(q);
              const auto& point_vectorized = phi.quadrature_point(q);
              auto u_m                     = u;
              auto grad_u_m                = grad_u;
              for(unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
                Point<dim> point;
                for(unsigned int d = 0; d < dim; ++d)
                  point[d] = point_vectorized[d][v];

                u_m[1][v] = -u_m[1][v];

                grad_u_m[0][0][v] = -grad_u_m[0][0][v];
                grad_u_m[0][1][v] = -grad_u_m[0][1][v];
              }

              const auto& visc = viscosity.value(phi_y_plus.get_value(q), grad_u, dx, Re);

              phi.submit_value(a33*(-(0.5*(visc*grad_u + viscosity.value(phi_y_plus.get_value(q), grad_u_m, dx, Re)*grad_u_m))*n_plus +
                              visc*coef_jump*(u - u_m)) + a33*outer_product(0.5*(u + u_m), phi_extr.get_value(q))*n_plus +
                              a33*0.5*lambda*(u - u_m), q);
              phi.submit_gradient(-theta_v*a33*visc*(outer_product(u - u_m, n_plus) + outer_product(n_plus, u - u_m)), q);
            }
            phi.integrate(true, true);
            diagonal[i] = phi.get_dof_value(i);
          }
          for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
            phi.submit_dof_value(diagonal[i], i);
          phi.distribute_local_to_global(dst);
        }
      }
    }
  }


  // Now we consider the pressure related bilinear forms. We first assemble diagonal cell term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_diagonal_cell_term_pressure(const MatrixFree<dim, Number>&               data,
                                       Vec&                                         dst,
                                       const unsigned int&                          ,
                                       const std::pair<unsigned int, unsigned int>& cell_range) const {
    FEEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number> phi(data, 1, 1);

    AlignedVector<VectorizedArray<Number>> diagonal(phi.dofs_per_component); /*--- Here we are using dofs_per_component but
                                                                                   it coincides with dofs_per_cell since it is
                                                                                   scalar finite element space ---*/

    const double coeff = (TR_BDF2_stage == 1) ? 1e6*gamma*dt*gamma*dt : 1e6*(1.0 - gamma)*dt*(1.0 - gamma)*dt;

    /*--- Loop over all cells in the range ---*/
    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
      phi.reinit(cell);

      /*--- Loop over all dofs ---*/
      for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
        for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
          phi.submit_dof_value(VectorizedArray<Number>(), j); /*--- We set all dofs to zero ---*/
        phi.submit_dof_value(make_vectorized_array<Number>(1.0), i); /*--- Now we set the current one to 1; since it is scalar,
                                                                           we can directly use 'make_vectorized_array' without
                                                                           relying on 'Tensor' ---*/
        phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

        /*--- Loop over quadrature points ---*/
        for(unsigned int q = 0; q < phi.n_q_points; ++q) {
          phi.submit_value(1.0/coeff*phi.get_value(q), q);
          phi.submit_gradient(phi.get_gradient(q), q);
        }
        phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        diagonal[i] = phi.get_dof_value(i);
      }
      for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
        phi.submit_dof_value(diagonal[i], i);

      phi.distribute_local_to_global(dst);
    }
  }


  // The following function assembles diagonal face term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_diagonal_face_term_pressure(const MatrixFree<dim, Number>&               data,
                                       Vec&                                         dst,
                                       const unsigned int&                          ,
                                       const std::pair<unsigned int, unsigned int>& face_range) const {
    FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number> phi_p(data, true, 1, 1),
                                                                   phi_m(data, false, 1, 1);

    AssertDimension(phi_p.dofs_per_component, phi_m.dofs_per_component);
    AlignedVector<VectorizedArray<Number>> diagonal_p(phi_p.dofs_per_component),
                                           diagonal_m(phi_m.dofs_per_component); /*--- Again, we just assert for safety that dimension
                                                                                       match, in the sense that we have selected
                                                                                       the proper space ---*/

    /*--- Loop over all faces ---*/
    for(unsigned int face = face_range.first; face < face_range.second; ++face) {
      phi_p.reinit(face);
      phi_m.reinit(face);

      const auto coef_jump = C_p*0.5*(std::abs((phi_p.get_normal_vector(0)*phi_p.inverse_jacobian(0))[dim - 1]) +
                                      std::abs((phi_m.get_normal_vector(0)*phi_m.inverse_jacobian(0))[dim - 1]));

      /*--- Loop over all dofs ---*/
      for(unsigned int i = 0; i < phi_p.dofs_per_component; ++i) {
        for(unsigned int j = 0; j < phi_p.dofs_per_component; ++j) {
          phi_p.submit_dof_value(VectorizedArray<Number>(), j);
          phi_m.submit_dof_value(VectorizedArray<Number>(), j);
        }
        phi_p.submit_dof_value(make_vectorized_array<Number>(1.0), i);
        phi_m.submit_dof_value(make_vectorized_array<Number>(1.0), i);
        phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

        /*--- Loop over all quadrature points to compute the integral ---*/
        for(unsigned int q = 0; q < phi_p.n_q_points; ++q) {
          const auto& n_plus        = phi_p.get_normal_vector(q);

          const auto& avg_grad_pres = 0.5*(phi_p.get_gradient(q) + phi_m.get_gradient(q));
          const auto& jump_pres     = phi_p.get_value(q) - phi_m.get_value(q);

          phi_p.submit_value(-scalar_product(avg_grad_pres, n_plus) + coef_jump*jump_pres, q);
          phi_m.submit_value(scalar_product(avg_grad_pres, n_plus) - coef_jump*jump_pres, q);
          phi_p.submit_gradient(-theta_p*0.5*jump_pres*n_plus, q);
          phi_m.submit_gradient(-theta_p*0.5*jump_pres*n_plus, q);
        }
        phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        diagonal_p[i] = phi_p.get_dof_value(i);
        phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        diagonal_m[i] = phi_m.get_dof_value(i);
      }
      for(unsigned int i = 0; i < phi_p.dofs_per_component; ++i) {
        phi_p.submit_dof_value(diagonal_p[i], i);
        phi_m.submit_dof_value(diagonal_m[i], i);
      }
      phi_p.distribute_local_to_global(dst);
      phi_m.distribute_local_to_global(dst);
    }
  }


  // Eventually, we assemble diagonal boundary term for the pressure
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  assemble_diagonal_boundary_term_pressure(const MatrixFree<dim, Number>&               data,
                                           Vec&                                         dst,
                                           const unsigned int&                          ,
                                           const std::pair<unsigned int, unsigned int>& face_range) const {
    FEFaceEvaluation<dim, fe_degree_p, n_q_points_1d_p, 1, Number> phi(data, true, 1, 1);

    AlignedVector<VectorizedArray<Number>> diagonal(phi.dofs_per_component);

    for(unsigned int face = face_range.first; face < face_range.second; ++face) {
      const auto boundary_id = data.get_boundary_id(face);

      if(boundary_id == 1) {
        phi.reinit(face);

        const auto coef_jump = C_p*std::abs((phi.get_normal_vector(0)*phi.inverse_jacobian(0))[dim - 1]);

        for(unsigned int i = 0; i < phi.dofs_per_component; ++i) {
          for(unsigned int j = 0; j < phi.dofs_per_component; ++j)
            phi.submit_dof_value(VectorizedArray<Number>(), j);
          phi.submit_dof_value(make_vectorized_array<Number>(1.0), i);
          phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

          for(unsigned int q = 0; q < phi.n_q_points; ++q) {
            const auto& n_plus    = phi.get_normal_vector(q);

            const auto& grad_pres = phi.get_gradient(q);
            const auto& pres      = phi.get_value(q);

            phi.submit_value(-scalar_product(grad_pres, n_plus) + 2.0*coef_jump*pres , q);
            phi.submit_normal_derivative(-theta_p*pres, q);
          }
          phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
          diagonal[i] = phi.get_dof_value(i);
        }
        for(unsigned int i = 0; i < phi.dofs_per_component; ++i)
          phi.submit_dof_value(diagonal[i], i);
        phi.distribute_local_to_global(dst);
      }
    }
  }


  // Put together all previous steps. We create a dummy auxliary vector that serves for the src input argument in
  // the previous functions that as we have seen before is unused. Then everything is done by the 'loop' function
  // and it is saved in the field 'inverse_diagonal_entries' already present in the base class. Anyway since there is
  // only one field, we need to resize properly depending on whether we are considering the velocity or the pressure.
  //
  template<int dim, int fe_degree_p, int fe_degree_v, int n_q_points_1d_p, int n_q_points_1d_v, typename Vec>
  void NavierStokesProjectionOperator<dim, fe_degree_p, fe_degree_v, n_q_points_1d_p, n_q_points_1d_v, Vec>::
  compute_diagonal() {
    Assert(NS_stage == 1 || NS_stage == 2, ExcInternalError());

    this->inverse_diagonal_entries.reset(new DiagonalMatrix<Vec>());
    auto& inverse_diagonal = this->inverse_diagonal_entries->get_vector();

    if(NS_stage == 1) {
      ::MatrixFreeTools::compute_diagonal<dim, Number, VectorizedArray<Number>>
      (*(this->data),
       inverse_diagonal,
       [&](const auto& data, auto& dst, const auto& src, const auto& cell_range) {
         (this->assemble_diagonal_cell_term_velocity)(data, dst, src, cell_range);
       },
       [&](const auto& data, auto& dst, const auto& src, const auto& face_range) {
         (this->assemble_diagonal_face_term_velocity)(data, dst, src, face_range);
       },
       [&](const auto& data, auto& dst, const auto& src, const auto& boundary_range) {
         (this->assemble_diagonal_boundary_term_velocity)(data, dst, src, boundary_range);
       },
       0);
    }
    else if(NS_stage == 2) {
        ::MatrixFreeTools::compute_diagonal<dim, Number, VectorizedArray<Number>>
      (*(this->data),
       inverse_diagonal,
       [&](const auto& data, auto& dst, const auto& src, const auto& cell_range) {
         (this->assemble_diagonal_cell_term_pressure)(data, dst, src, cell_range);
       },
       [&](const auto& data, auto& dst, const auto& src, const auto& face_range) {
         (this->assemble_diagonal_face_term_pressure)(data, dst, src, face_range);
       },
       [&](const auto& data, auto& dst, const auto& src, const auto& boundary_range) {
         (this->assemble_diagonal_boundary_term_pressure)(data, dst, src, boundary_range);
       },
       1);
    }

    for(unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i) {
      Assert(inverse_diagonal.local_element(i) != 0.0,
             ExcMessage("No diagonal entry in a definite operator should be zero"));
      inverse_diagonal.local_element(i) = 1.0/inverse_diagonal.local_element(i);
    }
  }


  // @sect{The <code>NavierStokesProjection</code> class}

  // Now we are ready for the main class of the program. It implements the calls to the various steps
  // of the projection method for Navier-Stokes equations.
  //
  template<int dim>
  class NavierStokesProjection {
  public:
    NavierStokesProjection(RunTimeParameters::Data_Storage& data);

    void run(const bool verbose = false, const unsigned int output_interval = 10);

  protected:
    const double t_0;
    const double T;
    const double gamma;         //--- TR-BDF2 parameter
    unsigned int TR_BDF2_stage; //--- Flag to check at which current stage of TR-BDF2 are
    const double Re;
    const double Cs2;
    double       dt;

    EquationData::Velocity<dim> vel_init;
    EquationData::Pressure<dim> pres_init; /*--- Instance of 'Velocity' and 'Pressure' classes to initialize. ---*/

    parallel::distributed::Triangulation<dim> triangulation;

    /*--- parameters to find nearest boundary vertex (for calculation of y+) ---*/
    Triangulation<dim> serial_triangulation;
    std::map<typename Triangulation<dim>::active_cell_iterator, int> cell_to_nearest_boundary_point;
    std::vector<std::vector<int>> cell_to_nearest_boundary_point_all;
    std::vector<std::vector<bool>> process_owns_boundary_point_all;
    std::map<int, int> owned_boundary_points;
    std::map<int, int> boundary_points_to_rank;

    /*--- Finite Element spaces ---*/
    FESystem<dim> fe_velocity;
    FESystem<dim> fe_pressure;
    FESystem<dim> fe_deltas;

    /*--- Handler for dofs ---*/
    DoFHandler<dim> dof_handler_velocity;
    DoFHandler<dim> dof_handler_pressure;
    DoFHandler<dim> dof_handler_deltas;

    /*--- Quadrature formulas for velocity and pressure, respectively ---*/
    QGauss<dim> quadrature_pressure;
    QGauss<dim> quadrature_velocity;


    /*--- Now we define all the vectors for the solution. We start from the pressure
          with p^n, p^(n+gamma) and a vector for rhs ---*/
    LinearAlgebra::distributed::Vector<double> pres_n;
    LinearAlgebra::distributed::Vector<double> pres_int;
    LinearAlgebra::distributed::Vector<double> rhs_p;

    /*--- Next, we move to the velocity, with u^n, u^(n-1), u^(n+gamma/2),
          u^(n+gamma) and other two auxiliary vectors as well as the rhs ---*/
    LinearAlgebra::distributed::Vector<double> u_n;
    LinearAlgebra::distributed::Vector<double> u_n_minus_1;
    LinearAlgebra::distributed::Vector<double> u_extr;
    LinearAlgebra::distributed::Vector<double> u_n_gamma;
    LinearAlgebra::distributed::Vector<double> u_star;
    LinearAlgebra::distributed::Vector<double> u_tmp;
    LinearAlgebra::distributed::Vector<double> rhs_u;
    LinearAlgebra::distributed::Vector<double> grad_pres_int;

    LinearAlgebra::distributed::Vector<double> deltas;
    LinearAlgebra::distributed::Vector<double> y_plus;
    LinearAlgebra::distributed::Vector<double> boundary_distance;

    LinearAlgebra::distributed::Vector<double> artificial_force;

    /*--- Variables for statistics ---*/
    std::vector<Point<dim>> obstacle_points;
    std::vector<Point<dim>> horizontal_wake_points;
    std::vector<Point<dim>> vertical_profile_points1;
    std::vector<Point<dim>> vertical_profile_points2;
    std::vector<Point<dim>> vertical_profile_points3;

    std::vector<double> avg_pressure;
    std::vector<double> avg_stress;

    std::vector<Vector<double>> avg_horizontal_velocity;
    std::vector<Vector<double>> avg_vertical_velocity1;
    std::vector<Vector<double>> avg_vertical_velocity2;
    std::vector<Vector<double>> avg_vertical_velocity3;

    Vector<double> Linfty_error_per_cell_vel;

    DeclException2(ExcInvalidTimeStep,
                   double,
                   double,
                   << " The time step " << arg1 << " is out of range."
                   << std::endl
                   << " The permitted range is (0," << arg2 << "]");

    void import_triangulation(const unsigned int n_refines, std::string filename, double x_start, double y_start );

    void create_triangulation_2D(const unsigned int n_refines);

    void create_triangulation_3D(const unsigned int n_refines);

    void setup_dofs();

    void initialize();

    void interpolate_velocity();

    void diffusion_step();

    void projection_step();

    void project_grad(const unsigned int flag);

    double get_maximal_velocity();

    double get_maximal_difference_velocity();

    void output_results(const unsigned int step);

    void output_statistics(Point<dim> center);

    void refine_mesh();

    void interpolate_max_res(const unsigned int level);

    void save_max_res();

  private:
    void compute_lift_and_drag();

    void initialize_points_around_obstacle(const unsigned int n_points, Point<dim> start, double dx);

    std::vector<Point<dim>> initialize_profile_points(double angle, double spacing, Point<dim> start_point,  Point<dim> end_point);

    void initialize_nearest_boundary_point_mapping();

    void compute_pressure_avg_over_boundary(int n, double height = 0.0, int n_points = 1);

    void compute_stress_avg_over_boundary(int n, Point<dim> center, double object_length, double lower_boundary, double upper_boundary);

    void compute_lipschitz_number();

    void compute_velocity_avg(int n, std::vector<Point<dim>>& points, std::vector<Vector<double>>& avg_velocity);

    void compute_artificial_force(LinearAlgebra::distributed::Vector<double> & vel);

    void compute_y_plus(LinearAlgebra::distributed::Vector<double> & vel, double lower_boundary, double upper_boundary, Point<dim> center, double object_length);

    /*--- Technical member to handle the various steps ---*/
    std::shared_ptr<MatrixFree<dim, double>> matrix_free_storage;

    /*--- Now we need an instance of the class implemented before with the weak form ---*/
    NavierStokesProjectionOperator<dim, EquationData::degree_p, EquationData::degree_p + 1,
                                   EquationData::degree_p + 1, EquationData::degree_p + 2,
                                   LinearAlgebra::distributed::Vector<double>> navier_stokes_matrix;

    /*--- This is an instance for geometric multigrid preconditioner ---*/
    MGLevelObject<NavierStokesProjectionOperator<dim, EquationData::degree_p, EquationData::degree_p + 1,
                                                 EquationData::degree_p + 1, EquationData::degree_p + 2,
                                                 LinearAlgebra::distributed::Vector<float>>> mg_matrices;

    /*--- Here we define two 'AffineConstraints' instance, one for each finite element space.
          This is just a technical issue, due to MatrixFree requirements. In general
          this class is used to impose boundary conditions (or any kind of constraints), but in this case, since
          we are using a weak imposition of bcs, everything is already in the weak forms and so these instances
          will be default constructed ---*/
    AffineConstraints<double> constraints_velocity,
                              constraints_pressure,
                              constraints_deltas;

    /*--- Now a bunch of variables handled by 'ParamHandler' introduced at the beginning of the code ---*/
    unsigned int max_its;
    double       eps;
    double       tolerance_fixed_point;

    unsigned int n_refines;
    unsigned int max_loc_refinements;
    unsigned int min_loc_refinements;
    unsigned int refinement_iterations;
    bool         import_mesh;
    bool         no_slip;

    std::string  saving_dir;

    bool         restart,
                 save_for_restart;
    unsigned int step_restart;
    double       time_restart;
    bool         as_initial_conditions;
    bool         modify_Reynolds;

    /*--- Finally, some output related streams ---*/
    ConditionalOStream pcout;

    std::ofstream      time_out;
    ConditionalOStream ptime_out;
    TimerOutput        time_table;

    std::ofstream output_n_dofs_velocity;
    std::ofstream output_n_dofs_pressure;

    std::ofstream output_lift;
    std::ofstream output_drag;

    std::ofstream output_avg_pressure;
    std::ofstream output_avg_stress;
    std::ofstream output_Cp;
    std::ofstream output_Cf;
    std::ofstream output_lipschitz;
    std::ofstream out_vel_hor;
    std::ofstream out_vel_ver1;
    std::ofstream out_vel_ver2;
    std::ofstream out_vel_ver3;
  };


  // In the constructor, we just read all the data from the
  // <code>Data_Storage</code> object that is passed as an argument, verify that
  // the data we read are reasonable and, finally, create the triangulation and
  // load the initial data.
  //
  template<int dim>
  NavierStokesProjection<dim>::NavierStokesProjection(RunTimeParameters::Data_Storage& data):
    t_0(data.initial_time),
    T(data.final_time),
    gamma(2.0 - std::sqrt(2.0)),  //--- Save also in the NavierStokes class the TR-BDF2 parameter value
    TR_BDF2_stage(1),             //--- Initialize the flag for the TR_BDF2 stage
    Re(data.Reynolds),
    Cs2(data.Cs2),
    dt(data.dt),
    vel_init(data.initial_time),
    pres_init(data.initial_time),
    triangulation(MPI_COMM_WORLD, parallel::distributed::Triangulation<dim>::limit_level_difference_at_vertices,
                  parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy),
    fe_velocity(FE_DGQ<dim>(EquationData::degree_p + 1), dim),
    fe_pressure(FE_DGQ<dim>(EquationData::degree_p), 1),
    fe_deltas(FE_DGQ<dim>(0), 1),
    dof_handler_velocity(triangulation),
    dof_handler_pressure(triangulation),
    dof_handler_deltas(triangulation),
    quadrature_pressure(EquationData::degree_p + 1),
    quadrature_velocity(EquationData::degree_p + 2),
    navier_stokes_matrix(data),
    max_its(data.max_iterations),
    eps(data.eps),
    tolerance_fixed_point(data.tolerance_fixed_point),
    n_refines(data.n_refines),
    import_mesh(data.import_mesh),
    no_slip(data.no_slip),
    max_loc_refinements(data.max_loc_refinements),
    min_loc_refinements(data.min_loc_refinements),
    refinement_iterations(data.refinement_iterations),
    saving_dir(data.dir),
    restart(data.restart),
    save_for_restart(data.save_for_restart),
    as_initial_conditions(data.as_initial_conditions),
    step_restart(data.step_restart),
    time_restart(data.time_restart),
    modify_Reynolds(data.modify_Reynolds),
    pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0),
    time_out("./" + data.dir + "/time_analysis_" +
             Utilities::int_to_string(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)) + "proc.dat"),
    ptime_out(time_out, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0),
    time_table(ptime_out, TimerOutput::summary, TimerOutput::cpu_and_wall_times),
    output_n_dofs_velocity("./" + data.dir + "/n_dofs_velocity.dat", std::ofstream::out),
    output_n_dofs_pressure("./" + data.dir + "/n_dofs_pressure.dat", std::ofstream::out),
    output_lift("./" + data.dir + "/lift.dat", std::ofstream::out),
    output_drag("./" + data.dir + "/drag.dat", std::ofstream::out),
    output_lipschitz("./" + data.dir + "/lipschitz.dat", std::ofstream::out)  {

      if(EquationData::degree_p < 1) {
        pcout
        << " WARNING: The chosen pair of finite element spaces is not stable."
        << std::endl
        << " The obtained results will be nonsense" << std::endl;
      }

      AssertThrow(!((dt <= 0.0) || (dt > 0.5*T)), ExcInvalidTimeStep(dt, 0.5*T));

      matrix_free_storage = std::make_shared<MatrixFree<dim, double>>();

      if(import_mesh){
        import_triangulation(n_refines, "unstr_sqcyl_coarse.msh", -10.0, -10.0);
      }
      else{
        if(dim == 2){
          create_triangulation_2D(n_refines);
        }
        else if(dim == 3){
          create_triangulation_3D(n_refines);
        }
      }

      setup_dofs();
      initialize();
  }

  template<int dim>
  void NavierStokesProjection<dim>::import_triangulation(const unsigned int n_refines, std::string filename, double x_start, double y_start){
    TimerOutput::Scope t(time_table, "Import triangulation");
    triangulation.clear();
    serial_triangulation.clear();

    /*--- parallel distributed triangulation ---*/
    GridIn<dim> gridin;
    gridin.attach_triangulation(triangulation);
    std::ifstream f(filename);
    gridin.read_msh(f);

    // /*--- Set boundary IDs ---*/
    // triangulation.set_all_manifold_ids_on_boundary(0, 0);
    // triangulation.set_all_manifold_ids_on_boundary(1, 1);
    // triangulation.set_all_manifold_ids_on_boundary(2, 2);
    // triangulation.set_all_manifold_ids_on_boundary(3, 3);

    /*--- Set boundary id ---*/
    for(const auto& face : triangulation.active_face_iterators()) {
      if(face->at_boundary()) {
        const Point<dim> center = face->center();
        // left side
        if(std::abs(center[0] - x_start) < 1e-10)
          face->set_boundary_id(0);
        // right side
        else if(std::abs(center[0] - (30.0+x_start)) < 1e-10)
          face->set_boundary_id(1);
        // cylinder boundary
        else if(center[0] < x_start + 10.5 + 1e-10 && center[0] > x_start + 9.5 - 1e-10 &&
                  center[1] < y_start + 10.5 + 1e-10 && center[1] > y_start + 9.5 - 1e-10)
          face->set_boundary_id(2);
        // sides of channel
        else {
          Assert(std::abs(center[1] - y_start) < 1.0e-10 ||
                std::abs(center[1] - (20.0+y_start)) < 1.0e-10,
                ExcInternalError());
          face->set_boundary_id(3);
        }
      }
    }

    /*--- We strongly advice to check the documentation to verify the meaning of all input parameters. ---*/
    if(restart) {
      triangulation.load("./" + saving_dir + "/solution_ser-" + Utilities::int_to_string(step_restart, 5));
    }
    else {
      pcout << "Number of refines = " << n_refines << std::endl;
      triangulation.refine_global(n_refines);
    }

    /*--- single triangulation ---*/
    GridIn<dim> sgridin;
    sgridin.attach_triangulation(serial_triangulation);
    std::ifstream sf(filename);
    sgridin.read_msh(sf);

    // /*--- Set boundary IDs ---*/
    // serial_triangulation.set_all_manifold_ids_on_boundary(0, 0);
    // serial_triangulation.set_all_manifold_ids_on_boundary(1, 1);
    // serial_triangulation.set_all_manifold_ids_on_boundary(2, 2);
    // serial_triangulation.set_all_manifold_ids_on_boundary(3, 3);

    /*--- Set boundary id ---*/
    for(const auto& face : serial_triangulation.active_face_iterators()) {
      if(face->at_boundary()) {
        const Point<dim> center = face->center();
        // left side
        if(std::abs(center[0] - x_start) < 1e-10)
          face->set_boundary_id(0);
        // right side
        else if(std::abs(center[0] - (30.0+x_start)) < 1e-10)
          face->set_boundary_id(1);
        // cylinder boundary
        else if(center[0] < x_start + 10.5 + 1e-10 && center[0] > x_start + 9.5 - 1e-10 &&
                  center[1] < y_start + 10.5 + 1e-10 && center[1] > y_start + 9.5 - 1e-10)
          face->set_boundary_id(2);
        // sides of channel
        else {
          Assert(std::abs(center[1] - y_start) < 1.0e-10 ||
                std::abs(center[1] - (20.0+y_start)) < 1.0e-10,
                ExcInternalError());
          face->set_boundary_id(3);
        }
      }
    }

    serial_triangulation.refine_global(n_refines);
  }

  // The method that creates the triangulation and refines it the needed number
  // of times.
  //
  template<int dim>
  void NavierStokesProjection<dim>::create_triangulation_2D(const unsigned int n_refines) {
    TimerOutput::Scope t(time_table, "Create triangulation");

    triangulation.clear();
    serial_triangulation.clear();

    /*--- parallel distributed triangulation ---*/
    parallel::distributed::Triangulation<dim> tria1(MPI_COMM_WORLD),
                                              tria2(MPI_COMM_WORLD),
                                              tria3(MPI_COMM_WORLD),
                                              tria4(MPI_COMM_WORLD),
                                              tria5(MPI_COMM_WORLD),
                                              tria6(MPI_COMM_WORLD),
                                              tria7(MPI_COMM_WORLD),
                                              tria8(MPI_COMM_WORLD),
                                              tria9(MPI_COMM_WORLD),
                                              tria10(MPI_COMM_WORLD),
                                              tria11(MPI_COMM_WORLD),
                                              tria12(MPI_COMM_WORLD);

    GridGenerator::subdivided_hyper_rectangle(tria1, {10, 16},
                                              Point<dim>(0.0, 10.7),
                                              Point<dim>(9.3, 20.0));
    GridGenerator::subdivided_hyper_rectangle(tria2, {14, 16},
                                              Point<dim>(9.3, 10.7),
                                              Point<dim>(10.7, 20.0));
    GridGenerator::subdivided_hyper_rectangle(tria3, {20, 16},
                                              Point<dim>(10.7, 10.7),
                                              Point<dim>(30.0, 20.0));
    GridGenerator::subdivided_hyper_rectangle(tria4, {10, 14},
                                              Point<dim>(0.0, 9.3),
                                              Point<dim>(9.3, 10.7));
    GridGenerator::subdivided_hyper_rectangle(tria5, {14, 2},
                                              Point<dim>(9.3, 10.5),
                                              Point<dim>(10.7, 10.7));
    GridGenerator::subdivided_hyper_rectangle(tria6, {2, 10},
                                              Point<dim>(9.3, 9.5),
                                              Point<dim>(9.5, 10.5));
    GridGenerator::subdivided_hyper_rectangle(tria7, {2, 10},
                                              Point<dim>(10.5, 9.5),
                                              Point<dim>(10.7, 10.5));
    GridGenerator::subdivided_hyper_rectangle(tria8, {14, 2},
                                              Point<dim>(9.3, 9.3),
                                              Point<dim>(10.7, 9.5));
    GridGenerator::subdivided_hyper_rectangle(tria9, {20, 14},
                                              Point<dim>(10.7, 9.3),
                                              Point<dim>(30.0, 10.7));
    GridGenerator::subdivided_hyper_rectangle(tria10, {10, 16},
                                              Point<dim>(0.0, 0.0),
                                              Point<dim>(9.3, 9.3));
    GridGenerator::subdivided_hyper_rectangle(tria11, {14, 16},
                                              Point<dim>(9.3, 0.0),
                                              Point<dim>(10.7, 9.3));
    GridGenerator::subdivided_hyper_rectangle(tria12, {20, 16},
                                              Point<dim>(10.7, 0.0),
                                              Point<dim>(30.0, 9.3));
    GridGenerator::merge_triangulations({&tria1, &tria2, &tria3, &tria4, &tria5, &tria6, &tria7, &tria8, &tria9, &tria10, &tria11, &tria12},
                                         triangulation, 1e-8, true);

    /*--- Set boundary id ---*/
    for(const auto& face : triangulation.active_face_iterators()) {
      if(face->at_boundary()) {
        const Point<dim> center = face->center();
        // left side
        if(std::abs(center[0] - 0.0) < 1e-10)
          face->set_boundary_id(0);
        // right side
        else if(std::abs(center[0] - 30.0) < 1e-10)
          face->set_boundary_id(1);
        // cylinder boundary
        else if(center[0] < 10.5 + 1e-10 && center[0] > 9.5 - 1e-10 &&
                  center[1] < 10.5 + 1e-10 && center[1] > 9.5 - 1e-10)
          face->set_boundary_id(2);
        // sides of channel
        else {
          Assert(std::abs(center[1] - 0.00) < 1.0e-10 ||
                std::abs(center[1] - 20.0) < 1.0e-10,
                ExcInternalError());
          face->set_boundary_id(3);
        }
      }
    }
    
    /*--- We strongly advice to check the documentation to verify the meaning of all input parameters. ---*/
    if(restart) {
      triangulation.load("./" + saving_dir + "/solution_ser-" + Utilities::int_to_string(step_restart, 5));
    }
    else {
      pcout << "Number of refines = " << n_refines << std::endl;
      triangulation.refine_global(n_refines);
    }

    /*--- single triangulation ---*/
    Triangulation<dim> stria1, stria2, stria3, stria4,
                       stria5, stria6, stria7, stria8,
                       stria9, stria10, stria11, stria12;

    GridGenerator::subdivided_hyper_rectangle(stria1, {10, 16},
                                              Point<dim>(0.0, 10.7),
                                              Point<dim>(9.3, 20.0));
    GridGenerator::subdivided_hyper_rectangle(stria2, {14, 16},
                                              Point<dim>(9.3, 10.7),
                                              Point<dim>(10.7, 20.0));
    GridGenerator::subdivided_hyper_rectangle(stria3, {20, 16},
                                              Point<dim>(10.7, 10.7),
                                              Point<dim>(30.0, 20.0));
    GridGenerator::subdivided_hyper_rectangle(stria4, {10, 14},
                                              Point<dim>(0.0, 9.3),
                                              Point<dim>(9.3, 10.7));
    GridGenerator::subdivided_hyper_rectangle(stria5, {14, 2},
                                              Point<dim>(9.3, 10.5),
                                              Point<dim>(10.7, 10.7));
    GridGenerator::subdivided_hyper_rectangle(stria6, {2, 10},
                                              Point<dim>(9.3, 9.5),
                                              Point<dim>(9.5, 10.5));
    GridGenerator::subdivided_hyper_rectangle(stria7, {2, 10},
                                              Point<dim>(10.5, 9.5),
                                              Point<dim>(10.7, 10.5));
    GridGenerator::subdivided_hyper_rectangle(stria8, {14, 2},
                                              Point<dim>(9.3, 9.3),
                                              Point<dim>(10.7, 9.5));
    GridGenerator::subdivided_hyper_rectangle(stria9, {20, 14},
                                              Point<dim>(10.7, 9.3),
                                              Point<dim>(30.0, 10.7));
    GridGenerator::subdivided_hyper_rectangle(stria10, {10, 16},
                                              Point<dim>(0.0, 0.0),
                                              Point<dim>(9.3, 9.3));
    GridGenerator::subdivided_hyper_rectangle(stria11, {14, 16},
                                              Point<dim>(9.3, 0.0),
                                              Point<dim>(10.7, 9.3));
    GridGenerator::subdivided_hyper_rectangle(stria12, {20, 16},
                                              Point<dim>(10.7, 0.0),
                                              Point<dim>(30.0, 9.3));
    GridGenerator::merge_triangulations({&stria1, &stria2, &stria3, &stria4, &stria5, &stria6, &stria7, &stria8, &stria9, &stria10, &stria11, &stria12},
                                         serial_triangulation, 1e-8, true);


    /*--- Set boundary id ---*/
    for(const auto& face : serial_triangulation.active_face_iterators()) {
      if(face->at_boundary()) {
        const Point<dim> center = face->center();

        // left side
        if(std::abs(center[0] - 0.0) < 1e-10)
          face->set_boundary_id(0);
        // right side
        else if(std::abs(center[0] - 30.0) < 1e-10)
          face->set_boundary_id(1);
        // cylinder boundary
        else if(center[0] < 10.5 + 1e-10 && center[0] > 9.5 - 1e-10 &&
                  center[1] < 10.5 + 1e-10 && center[1] > 9.5 - 1e-10)
          face->set_boundary_id(2);
        // sides of channel
        else {
          Assert(std::abs(center[1] - 0.00) < 1.0e-10 ||
                std::abs(center[1] - 20.0) < 1.0e-10,
                ExcInternalError());
          face->set_boundary_id(3);
        }
      }
    }
    serial_triangulation.refine_global(n_refines);
  }

  template<int dim>
  void NavierStokesProjection<dim>::create_triangulation_3D(const unsigned int n_refines) {
    TimerOutput::Scope t(time_table, "Create triangulation");

    triangulation.clear();
    serial_triangulation.clear();

    /*--- parallel distributed triangulation ---*/
    parallel::distributed::Triangulation<dim> tria1(MPI_COMM_WORLD),
                                                tria2(MPI_COMM_WORLD),
                                                tria3(MPI_COMM_WORLD),
                                                tria4(MPI_COMM_WORLD);

    GridGenerator::subdivided_hyper_rectangle(tria1, {60, 19, 3},
                                              Point<dim>(0.0, 0.0, 0.0),
                                              Point<dim>(30.0, 9.5, numbers::PI));
    GridGenerator::subdivided_hyper_rectangle(tria2, {19, 2, 3},
                                              Point<dim>(0.0, 9.5, 0.0),
                                              Point<dim>(9.5, 10.5, numbers::PI));
    GridGenerator::subdivided_hyper_rectangle(tria3, {39, 2, 3},
                                              Point<dim>(10.5, 9.5, 0.0),
                                              Point<dim>(30.0, 10.5, numbers::PI));
    GridGenerator::subdivided_hyper_rectangle(tria4, {60, 19, 3},
                                              Point<dim>(0.0, 10.5, 0.0),
                                              Point<dim>(30.0, 20.0, numbers::PI));
    GridGenerator::merge_triangulations({&tria1, &tria2, &tria3, &tria4},
                                        triangulation, 1e-8, true);

    /*--- Set boundary id ---*/
    for(const auto& face : triangulation.active_face_iterators()) {
      if(face->at_boundary()) {
        const Point<dim> center = face->center();

        // left side
        if(std::abs(center[0] - 0.0) < 1e-10)
          face->set_boundary_id(0);
        // right side
        else if(std::abs(center[0] - 30.0) < 1e-10)
          face->set_boundary_id(1);
        // sides of channel
        else if(std::abs(center[1] - 0.0) < 1.0e-10 || std::abs(center[1] - 20.0) < 1.0e-10)
          face->set_boundary_id(3);
        else if(std::abs(center[2] - 0.0) < 1e-10)
          face->set_boundary_id(4);
        else {
          Assert(std::abs(center[2] - numbers::PI) < 1.0e-10,
                ExcInternalError());
          face->set_boundary_id(5);
        }
      }
    }

    /*--- We strongly advice to check the documentation to verify the meaning of all input parameters. ---*/
    if(restart) {
      triangulation.load("./" + saving_dir + "/solution_ser-" + Utilities::int_to_string(step_restart, 5));
    }
    else {
      pcout << "Number of refines = " << n_refines << std::endl;
      triangulation.refine_global(n_refines);
    }

    /*--- single triangulation ---*/
    Triangulation<dim> stria1, stria2, stria3, stria4;

    GridGenerator::subdivided_hyper_rectangle(stria1, {60, 19, 3},
                                              Point<dim>(0.0, 0.0, 0.0),
                                              Point<dim>(30.0, 9.5, numbers::PI));
    GridGenerator::subdivided_hyper_rectangle(stria2, {19, 2, 3},
                                              Point<dim>(0.0, 9.5, 0.0),
                                              Point<dim>(9.5, 10.5, numbers::PI));
    GridGenerator::subdivided_hyper_rectangle(stria3, {39, 2, 3},
                                              Point<dim>(10.5, 9.5, 0.0),
                                              Point<dim>(30.0, 10.5, numbers::PI));
    GridGenerator::subdivided_hyper_rectangle(stria4, {60, 19, 3},
                                              Point<dim>(0.0, 10.5, 0.0),
                                              Point<dim>(30.0, 20.0, numbers::PI));
    GridGenerator::merge_triangulations({&stria1, &stria2, &stria3, &stria4},
                                        serial_triangulation, 1e-8, true);

    /*--- Set boundary id ---*/
    for(const auto& face : serial_triangulation.active_face_iterators()) {
      if(face->at_boundary()) {
        const Point<dim> center = face->center();

        // left side
        if(std::abs(center[0] - 0.0) < 1e-10)
          face->set_boundary_id(0);
        // right side
        else if(std::abs(center[0] - 30.0) < 1e-10)
          face->set_boundary_id(1);
        // sides of channel
        else if(std::abs(center[1] - 0.0) < 1.0e-10 || std::abs(center[1] - 20.0) < 1.0e-10)
          face->set_boundary_id(3);
        else if(std::abs(center[2] - 0.0) < 1e-10)
          face->set_boundary_id(4);
        else {
          Assert(std::abs(center[2] - numbers::PI) < 1.0e-10,
                ExcInternalError());
          face->set_boundary_id(5);
        }
      }
    }
    serial_triangulation.refine_global(n_refines);
  }

  // After creating the triangulation, it creates the mesh dependent
  // data, i.e. it distributes degrees of freedom, and
  // initializes the vectors that we will use.
  //
  template<int dim>
  void NavierStokesProjection<dim>::setup_dofs() {
    pcout << "Number of active cells: " << triangulation.n_global_active_cells() << std::endl;
    pcout << "Number of levels: "       << triangulation.n_global_levels()       << std::endl;

    /*--- Distribute dofs and prepare for multigrid ---*/
    dof_handler_velocity.distribute_dofs(fe_velocity);
    dof_handler_pressure.distribute_dofs(fe_pressure);
    dof_handler_deltas.distribute_dofs(fe_deltas);

    pcout << "dim (X_h) = " << dof_handler_velocity.n_dofs()
          << std::endl
          << "dim (M_h) = " << dof_handler_pressure.n_dofs()
          << std::endl
          << "Re        = " << Re << std::endl
          << "Cs2       = " << Cs2 << std::endl
          << std::endl;

    if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
      output_n_dofs_velocity << dof_handler_velocity.n_dofs() << std::endl;
      output_n_dofs_pressure << dof_handler_pressure.n_dofs() << std::endl;
    }

    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags                = (update_gradients | update_JxW_values |
                                                           update_quadrature_points | update_values);
    additional_data.mapping_update_flags_inner_faces    = (update_gradients | update_JxW_values | update_quadrature_points |
                                                           update_normal_vectors | update_values);
    additional_data.mapping_update_flags_boundary_faces = (update_gradients | update_JxW_values | update_quadrature_points |
                                                           update_normal_vectors | update_values);
    additional_data.tasks_parallel_scheme               = MatrixFree<dim, double>::AdditionalData::none;

    std::vector<const DoFHandler<dim>*> dof_handlers; /*--- Vector of dof_handlers to feed the 'MatrixFree'. Here the order
                                                            counts and enters into the game as parameter of FEEvaluation and
                                                            FEFaceEvaluation in the previous class ---*/
    dof_handlers.push_back(&dof_handler_velocity);
    dof_handlers.push_back(&dof_handler_pressure);
    dof_handlers.push_back(&dof_handler_deltas);

    constraints_velocity.clear();
    constraints_velocity.close();
    constraints_pressure.clear();
    constraints_pressure.close();
    constraints_deltas.clear();
    constraints_deltas.close();
    std::vector<const AffineConstraints<double>*> constraints;
    constraints.push_back(&constraints_velocity);
    constraints.push_back(&constraints_pressure);
    constraints.push_back(&constraints_deltas);

    std::vector<QGauss<1>> quadratures; /*--- We cannot directly use 'quadrature_velocity' and 'quadrature_pressure',
                                              because the 'MatrixFree' structure wants a quadrature formula for 1D
                                              (this is way the template parameter of the previous class was called 'n_q_points_1d_p'
                                               and 'n_q_points_1d_v' and the reason of '1' as QGauss template parameter). ---*/
    quadratures.push_back(QGauss<1>(EquationData::degree_p + 2));
    quadratures.push_back(QGauss<1>(EquationData::degree_p + 1));
    quadratures.push_back(QGauss<1>(1));

    /*--- Initialize the matrix-free structure and size properly the vectors. Here again the
          second input argument of the 'initialize_dof_vector' method depends on the order of 'dof_handlers' ---*/
    matrix_free_storage->reinit(MappingQ1<dim>(), dof_handlers, constraints, quadratures, additional_data);
    matrix_free_storage->initialize_dof_vector(u_star, 0);
    matrix_free_storage->initialize_dof_vector(rhs_u, 0);
    matrix_free_storage->initialize_dof_vector(u_n, 0);
    matrix_free_storage->initialize_dof_vector(u_extr, 0);
    matrix_free_storage->initialize_dof_vector(u_n_minus_1, 0);
    matrix_free_storage->initialize_dof_vector(u_n_gamma, 0);
    matrix_free_storage->initialize_dof_vector(u_tmp, 0);
    matrix_free_storage->initialize_dof_vector(grad_pres_int, 0);
    matrix_free_storage->initialize_dof_vector(artificial_force, 0);

    matrix_free_storage->initialize_dof_vector(pres_int, 1);
    matrix_free_storage->initialize_dof_vector(pres_n, 1);
    matrix_free_storage->initialize_dof_vector(rhs_p, 1);

    matrix_free_storage->initialize_dof_vector(deltas, 2);
    matrix_free_storage->initialize_dof_vector(y_plus, 2);
    matrix_free_storage->initialize_dof_vector(boundary_distance, 2);


    // initialize delta calculation
    for(const auto& cell: dof_handler_deltas.active_cell_iterators()) {
      if(cell->is_locally_owned()) {
        std::vector<types::global_dof_index> dof_indices(fe_deltas.dofs_per_cell);
        cell->get_dof_indices(dof_indices);
        for(unsigned int idx = 0; idx < dof_indices.size(); ++idx) {
          deltas(dof_indices[idx]) = std::pow(cell->measure(), 1./dim)/(EquationData::degree_p + 1);
        }
      }
    }

    navier_stokes_matrix.set_deltas(deltas);

    /*--- Initialize the multigrid structure. We dedicate ad hoc 'dof_handlers_mg' and 'constraints_mg' because
          we use float as type. Moreover we can initialize already with the index of the finite element of the pressure;
          anyway we need by requirement to declare also structures for the velocity for coherence (basically because
          the index of finite element space has to be the same, so the pressure has to be the second).---*/
    mg_matrices.clear_elements();
    dof_handler_velocity.distribute_mg_dofs();
    dof_handler_pressure.distribute_mg_dofs();
    dof_handler_deltas.distribute_mg_dofs();

    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels - 1);
    for(unsigned int level = 0; level < nlevels; ++level) {
      typename MatrixFree<dim, float>::AdditionalData additional_data_mg;
      additional_data_mg.tasks_parallel_scheme               = MatrixFree<dim, float>::AdditionalData::none;
      additional_data_mg.mapping_update_flags                = (update_gradients | update_JxW_values);
      additional_data_mg.mapping_update_flags_inner_faces    = (update_gradients | update_JxW_values);
      additional_data_mg.mapping_update_flags_boundary_faces = (update_gradients | update_JxW_values);
      additional_data_mg.mg_level = level;

      std::vector<const DoFHandler<dim>*> dof_handlers_mg;
      dof_handlers_mg.push_back(&dof_handler_velocity);
      dof_handlers_mg.push_back(&dof_handler_pressure);
      dof_handlers_mg.push_back(&dof_handler_deltas);

      std::vector<const AffineConstraints<float>*> constraints_mg;
      AffineConstraints<float> constraints_velocity_mg;
      constraints_velocity_mg.clear();
      constraints_velocity_mg.close();
      constraints_mg.push_back(&constraints_velocity_mg);
      AffineConstraints<float> constraints_pressure_mg;
      constraints_pressure_mg.clear();
      constraints_pressure_mg.close();
      constraints_mg.push_back(&constraints_pressure_mg);
      AffineConstraints<float> constraints_deltas_mg;
      constraints_deltas_mg.clear();
      constraints_deltas_mg.close();
      constraints_mg.push_back(&constraints_deltas_mg);

      std::shared_ptr<MatrixFree<dim, float>> mg_mf_storage_level(new MatrixFree<dim, float>());
      mg_mf_storage_level->reinit(MappingQ1<dim>(),dof_handlers_mg, constraints_mg, quadratures, additional_data_mg);
      const std::vector<unsigned int> tmp = {1};
      mg_matrices[level].initialize(mg_mf_storage_level, tmp, tmp);
      mg_matrices[level].set_dt(dt);
      mg_matrices[level].set_NS_stage(2);
    }

    Linfty_error_per_cell_vel.reinit(triangulation.n_active_cells());
  }


  // This method loads the initial data. It simply uses the class <code>Pressure</code> instance for the pressure
  // and the class <code>Velocity</code> instance for the velocity.
  //
  template<int dim>
  void NavierStokesProjection<dim>::initialize() {
    TimerOutput::Scope t(time_table, "Initialize pressure and velocity");

    if(restart) {
      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
      solution_transfer_velocity(dof_handler_velocity);
      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
      solution_transfer_pressure(dof_handler_pressure);

      u_n.zero_out_ghost_values();
      u_n_minus_1.zero_out_ghost_values();
      pres_n.zero_out_ghost_values();

      std::vector<LinearAlgebra::distributed::Vector<double>*> velocities;
      velocities.push_back(&u_n);
      velocities.push_back(&u_n_minus_1);

      solution_transfer_velocity.deserialize(velocities);
      solution_transfer_pressure.deserialize(pres_n);

      if(n_refines - (triangulation.n_global_levels() - 1) > 0) {
        std::vector<const LinearAlgebra::distributed::Vector<double>*> velocities_tmp;
        velocities_tmp.push_back(&u_n);
        velocities_tmp.push_back(&u_n_minus_1);

        solution_transfer_velocity.prepare_for_coarsening_and_refinement(velocities_tmp);
        solution_transfer_pressure.prepare_for_coarsening_and_refinement(pres_n);

        triangulation.refine_global(n_refines - (triangulation.n_global_levels() - 1));

        setup_dofs();

        LinearAlgebra::distributed::Vector<double> transfer_velocity,
                                                   transfer_velocity_minus_1,
                                                   transfer_pressure;
        transfer_velocity.reinit(u_n);
        transfer_velocity.zero_out_ghost_values();
        transfer_velocity_minus_1.reinit(u_n_minus_1);
        transfer_velocity_minus_1.zero_out_ghost_values();
        transfer_pressure.reinit(pres_n);
        transfer_pressure.zero_out_ghost_values();

        std::vector<LinearAlgebra::distributed::Vector<double>*> transfer_velocities;
        transfer_velocities.push_back(&transfer_velocity);
        transfer_velocities.push_back(&transfer_velocity_minus_1);
        solution_transfer_velocity.interpolate(transfer_velocities);
        transfer_velocity.update_ghost_values();
        transfer_velocity_minus_1.update_ghost_values();
        solution_transfer_pressure.interpolate(transfer_pressure);
        transfer_pressure.update_ghost_values();

        u_n         = transfer_velocity;
        u_n_minus_1 = transfer_velocity_minus_1;
        pres_n      = transfer_pressure;
      }
    }
    else {
      VectorTools::interpolate(dof_handler_pressure, pres_init, pres_n);

      VectorTools::interpolate(dof_handler_velocity, vel_init, u_n_minus_1);
      VectorTools::interpolate(dof_handler_velocity, vel_init, u_n);
    }
  }


  // This function computes the extrapolated velocity to be used in the momentum predictor
  //
  template<int dim>
  void NavierStokesProjection<dim>::interpolate_velocity() {
    TimerOutput::Scope t(time_table, "Interpolate velocity");

    //--- TR-BDF2 first step
    if(TR_BDF2_stage == 1) {
      u_extr.equ(1.0 + gamma/(2.0*(1.0 - gamma)), u_n);
      u_tmp.equ(gamma/(2.0*(1.0 - gamma)), u_n_minus_1);
      u_extr -= u_tmp;
    }
    //--- TR-BDF2 second step
    else {
      u_extr.equ(1.0 + (1.0 - gamma)/gamma, u_n_gamma);
      u_tmp.equ((1.0 - gamma)/gamma, u_n);
      u_extr -= u_tmp;
    }
  }


  // We are finally ready to solve the diffusion step.
  //
  template<int dim>
  void NavierStokesProjection<dim>::diffusion_step() {
    TimerOutput::Scope t(time_table, "Diffusion step");

    /*--- We first speicify that we want to deal with velocity dof_handler (index 0, since it is the first one
          in the 'dof_handlers' vector) ---*/
    const std::vector<unsigned int> tmp = {0};
    navier_stokes_matrix.initialize(matrix_free_storage, tmp, tmp);

    /*--- Next, we specify at we are at stage 1, namely the diffusion step ---*/
    navier_stokes_matrix.set_NS_stage(1);

    /*--- Now, we compute the right-hand side and we set the convective velocity. The necessity of 'set_u_extr' is
          that this quantity is required in the bilinear forms and we can't use a vector of src like on the right-hand side,
          so it has to be available ---*/
    if(TR_BDF2_stage == 1) {
      navier_stokes_matrix.vmult_rhs_velocity(rhs_u, {u_n, u_extr, pres_n, deltas, artificial_force, y_plus});
      navier_stokes_matrix.set_u_extr(u_extr);
      u_star = u_extr;
    }
    else {
      navier_stokes_matrix.vmult_rhs_velocity(rhs_u, {u_n, u_n_gamma, pres_int, u_extr, deltas, artificial_force, y_plus});
      navier_stokes_matrix.set_u_extr(u_extr);
      u_star = u_extr;
    }

    /*--- Build the linear solver; in this case we specifiy the maximum number of iterations and residual ---*/
    SolverControl solver_control(max_its, eps*rhs_u.l2_norm());
    SolverGMRES<LinearAlgebra::distributed::Vector<double>> gmres(solver_control);

    /*--- Build a Jacobi preconditioner and solve ---*/
    PreconditionJacobi<NavierStokesProjectionOperator<dim,
                                                      EquationData::degree_p,
                                                      EquationData::degree_p + 1,
                                                      EquationData::degree_p + 1,
                                                      EquationData::degree_p + 2,
                                                      LinearAlgebra::distributed::Vector<double>>> preconditioner;
    navier_stokes_matrix.compute_diagonal();
    preconditioner.initialize(navier_stokes_matrix);

    gmres.solve(navier_stokes_matrix, u_star, rhs_u, preconditioner);
  }


  // Next, we solve the projection step.
  //
  template<int dim>
  void NavierStokesProjection<dim>::projection_step() {
    TimerOutput::Scope t(time_table, "Projection step pressure");

    /*--- We start in the same way of 'diffusion_step': we first reinitialize with the index of FE space,
          we specify that this is the second stage and we compute the right-hand side ---*/
    const std::vector<unsigned int> tmp = {1};
    navier_stokes_matrix.initialize(matrix_free_storage, tmp, tmp);

    navier_stokes_matrix.set_NS_stage(2);

    if(TR_BDF2_stage == 1)
      navier_stokes_matrix.vmult_rhs_pressure(rhs_p, {u_star, pres_n});
    else
      navier_stokes_matrix.vmult_rhs_pressure(rhs_p, {u_star, pres_int});

    /*--- Build the linear solver (Conjugate Gradient in this case) ---*/
    SolverControl solver_control(max_its, eps*rhs_p.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);

    /*--- Build the preconditioner (as in step-37) ---*/
    MGTransferMatrixFree<dim, float> mg_transfer;
    mg_transfer.build(dof_handler_pressure);

    using SmootherType = PreconditionChebyshev<NavierStokesProjectionOperator<dim,
                                                                              EquationData::degree_p,
                                                                              EquationData::degree_p + 1,
                                                                              EquationData::degree_p + 1,
                                                                              EquationData::degree_p + 2,
                                                                              LinearAlgebra::distributed::Vector<float>>,
                                               LinearAlgebra::distributed::Vector<float>>;
    mg::SmootherRelaxation<SmootherType, LinearAlgebra::distributed::Vector<float>> mg_smoother;
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
    smoother_data.resize(0, triangulation.n_global_levels() - 1);
    for(unsigned int level = 0; level < triangulation.n_global_levels(); ++level) {
      if(level > 0) {
        smoother_data[level].smoothing_range     = 15.0;
        smoother_data[level].degree              = 3;
        smoother_data[level].eig_cg_n_iterations = 10;
      }
      else {
        smoother_data[0].smoothing_range     = 2e-2;
        smoother_data[0].degree              = numbers::invalid_unsigned_int;
        smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
      }
      mg_matrices[level].compute_diagonal();
      smoother_data[level].preconditioner = mg_matrices[level].get_matrix_diagonal_inverse();
    }
    mg_smoother.initialize(mg_matrices, smoother_data);

    PreconditionIdentity                                identity;
    SolverCG<LinearAlgebra::distributed::Vector<float>> cg_mg(solver_control);
    MGCoarseGridIterativeSolver<LinearAlgebra::distributed::Vector<float>,
                                SolverCG<LinearAlgebra::distributed::Vector<float>>,
                                NavierStokesProjectionOperator<dim,
                                                               EquationData::degree_p,
                                                               EquationData::degree_p + 1,
                                                               EquationData::degree_p + 1,
                                                               EquationData::degree_p + 2,
                                                               LinearAlgebra::distributed::Vector<float>>,
                                PreconditionIdentity> mg_coarse(cg_mg, mg_matrices[0], identity);

    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(mg_matrices);

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg(mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);

    PreconditionMG<dim,
                   LinearAlgebra::distributed::Vector<float>,
                   MGTransferMatrixFree<dim, float>> preconditioner(dof_handler_pressure, mg, mg_transfer);

    /*--- Solve the linear system ---*/
    if(TR_BDF2_stage == 1) {
      pres_int = pres_n;
      cg.solve(navier_stokes_matrix, pres_int, rhs_p, preconditioner);
    }
    else {
      pres_n = pres_int;
      cg.solve(navier_stokes_matrix, pres_n, rhs_p, preconditioner);
    }
  }


  // This implements the projection step for the gradient of pressure
  //
  template<int dim>
  void NavierStokesProjection<dim>::project_grad(const unsigned int flag) {
    TimerOutput::Scope t(time_table, "Gradient of pressure projection");

    /*--- The input parameter flag is used just to specify where we want to save the result ---*/
    AssertIndexRange(flag, 3);
    Assert(flag > 0, ExcInternalError());

    /*--- We need to select the dof handler related to the velocity since the result lives there ---*/
    const std::vector<unsigned int> tmp = {0};
    navier_stokes_matrix.initialize(matrix_free_storage, tmp, tmp);

    if(flag == 1)
      navier_stokes_matrix.vmult_grad_p_projection(rhs_u, pres_n);
    else if(flag == 2)
      navier_stokes_matrix.vmult_grad_p_projection(rhs_u, pres_int);

    /*--- We conventionally decide that the this corresponds to third stage ---*/
    navier_stokes_matrix.set_NS_stage(3);

    /*--- Solve the system ---*/
    SolverControl solver_control(max_its, 1e-12*rhs_u.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);
    cg.solve(navier_stokes_matrix, u_tmp, rhs_u, PreconditionIdentity());
  }


  // The following function is used in determining the maximal velocity
  // in order to compute the Courant number.
  //
  template<int dim>
  double NavierStokesProjection<dim>::get_maximal_velocity() {
    return u_n.linfty_norm();
  }


  // The following function is used in determining the maximal nodal difference
  // between old and current velocity value in order to see if we have reched steady-state.
  //
  template<int dim>
  double NavierStokesProjection<dim>::get_maximal_difference_velocity() {
    u_tmp = u_n;
    u_tmp -= u_n_minus_1;

    return u_tmp.linfty_norm();
  }


  // This method plots the current solution. The main difficulty is that we want
  // to create a single output file that contains the data for all velocity
  // components and the pressure. On the other hand, velocities and the pressure
  // live on separate DoFHandler objects, so we need to pay attention when we use
  // 'add_data_vector' to select the proper space.
  //
  template<int dim>
  void NavierStokesProjection<dim>::output_results(const unsigned int step) {
    TimerOutput::Scope t(time_table, "Output results");

    DataOut<dim> data_out;

    std::vector<std::string> velocity_names(dim, "v");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation_velocity(dim, DataComponentInterpretation::component_is_part_of_vector);
    u_n.update_ghost_values();
    data_out.add_data_vector(dof_handler_velocity, u_n, velocity_names, component_interpretation_velocity);
    pres_n.update_ghost_values();
    data_out.add_data_vector(dof_handler_pressure, pres_n, "p", {DataComponentInterpretation::component_is_scalar});

    boundary_distance.update_ghost_values();
    data_out.add_data_vector(dof_handler_deltas, boundary_distance, "d", {DataComponentInterpretation::component_is_scalar});

    y_plus.update_ghost_values();
    data_out.add_data_vector(dof_handler_deltas, y_plus, "y_plus", {DataComponentInterpretation::component_is_scalar});

    std::vector<std::string> velocity_names_old(dim, "v_old");
    u_n_minus_1.update_ghost_values();
    data_out.add_data_vector(dof_handler_velocity, u_n_minus_1, velocity_names_old, component_interpretation_velocity);

    /*--- Here we rely on the postprocessor we have built ---*/
    PostprocessorVorticity<dim> postprocessor;
    data_out.add_data_vector(dof_handler_velocity, u_n, postprocessor);

    data_out.build_patches(MappingQ1<dim>(),
                           EquationData::degree_p + 1, DataOut<dim>::curved_inner_cells);

    const std::string output = "./" + saving_dir + "/solution-" + Utilities::int_to_string(step, 5) + ".vtu";
    data_out.write_vtu_in_parallel(output, MPI_COMM_WORLD);

    /*--- Serialization ---*/
    if(save_for_restart) {
      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
      solution_transfer_velocity(dof_handler_velocity);
      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
      solution_transfer_pressure(dof_handler_pressure);

      u_n.update_ghost_values();
      u_n_minus_1.update_ghost_values();
      pres_n.update_ghost_values();

      std::vector<const LinearAlgebra::distributed::Vector<double>*> velocities;
      velocities.push_back(&u_n);
      velocities.push_back(&u_n_minus_1);
      solution_transfer_velocity.prepare_for_serialization(velocities);
      solution_transfer_pressure.prepare_for_serialization(pres_n);

      triangulation.save("./" + saving_dir + "/solution_ser-" + Utilities::int_to_string(step, 5));
    }
  }

  template<int dim>
  void NavierStokesProjection<dim>::output_statistics(Point<dim> center) {

    const double p_inf = 30.0;
    const double U_inf = 1.0;

    output_avg_pressure.close();
    output_avg_stress.close();
    output_Cf.close();
    out_vel_hor.close();
    out_vel_ver1.close();
    out_vel_ver2.close();
    out_vel_ver3.close();
    output_Cp.close();

    for(unsigned int rank = 0; rank < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++rank) {

      if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == rank) {
        if(rank == 0){
          output_avg_pressure.open("./" + saving_dir + "/avg_p.dat", std::ofstream::out | std::ofstream::trunc);
          output_avg_stress.open("./" + saving_dir + "/avg_stress.dat", std::ofstream::out | std::ofstream::trunc);
          output_Cf.open("./" + saving_dir + "/Cf.dat", std::ofstream::out | std::ofstream::trunc);
          output_Cp.open("./" + saving_dir + "/Cp.dat", std::ofstream::out | std::ofstream::trunc);
          out_vel_hor.open("./" + saving_dir + "/out_vel_hor.dat", std::ofstream::out | std::ofstream::trunc);
          out_vel_ver1.open("./" + saving_dir + "/out_vel_ver1.dat", std::ofstream::out | std::ofstream::trunc);
          out_vel_ver2.open("./" + saving_dir + "/out_vel_ver2.dat", std::ofstream::out | std::ofstream::trunc);
          out_vel_ver3.open("./" + saving_dir + "/out_vel_ver3.dat", std::ofstream::out | std::ofstream::trunc);
        }
        else{
          output_avg_pressure.open("./" + saving_dir + "/avg_p.dat", std::ios_base::app);
          output_avg_stress.open("./" + saving_dir + "/avg_stress.dat", std::ios_base::app);
          output_Cf.open("./" + saving_dir + "/Cf.dat", std::ios_base::app);
          output_Cp.open("./" + saving_dir + "/Cp.dat", std::ios_base::app);
          out_vel_hor.open("./" + saving_dir + "/out_vel_hor.dat", std::ios_base::app);
          out_vel_ver1.open("./" + saving_dir + "/out_vel_ver1.dat", std::ios_base::app);
          out_vel_ver2.open("./" + saving_dir + "/out_vel_ver2.dat", std::ios_base::app);
          out_vel_ver3.open("./" + saving_dir + "/out_vel_ver3.dat", std::ios_base::app);
        }

        /*--- Output average along cylinder boundary ---*/
        for(unsigned int i = 0; i < obstacle_points.size(); i++) {

          /*--- Output pressure average ---*/
          output_avg_pressure << obstacle_points[i][0] << " " << obstacle_points[i][1] << " " << avg_pressure[i] << std::endl;

          /*--- Output Cf average ---*/
          output_Cf << obstacle_points[i][0] << " " << obstacle_points[i][1] << " " << avg_stress[i] * 2. / (Re * U_inf * U_inf) << std::endl;

          /*--- Output stress average ---*/
          output_avg_stress << obstacle_points[i][0] << " " << obstacle_points[i][1] << " " << avg_stress[i] << std::endl;

          /*--- Output Cp average ---*/
          output_Cp << obstacle_points[i][0] << " " << obstacle_points[i][1] << " " << 2.0 * (avg_pressure[i] - p_inf) / (U_inf*U_inf) << std::endl;
        }

        /*--- Output average velocity horizontal wake points ---*/
        for(unsigned int i = 0; i < horizontal_wake_points.size(); i++) {
          out_vel_hor  << horizontal_wake_points[i][0] << " "
                       << avg_horizontal_velocity[i][0] << " "
                       << avg_horizontal_velocity[i][1] << std::endl;
        }

        /*--- Output average velocity vertical points 1 ---*/
        for(unsigned int i = 0; i < vertical_profile_points1.size(); i++) {
          out_vel_ver1 << vertical_profile_points1[i][1] << " "
                       << avg_vertical_velocity1[i][0] << " "
                       << avg_vertical_velocity1[i][1] << std::endl;
        }

        /*--- Output average velocity vertical points 2 ---*/
        for(unsigned int i = 0; i < vertical_profile_points2.size(); i++) {
          out_vel_ver2 << vertical_profile_points2[i][1] << " "
                       << avg_vertical_velocity2[i][0] << " "
                       << avg_vertical_velocity2[i][1] << std::endl;
        }

        /*--- Output average velocity vertical points 3 ---*/
        for(unsigned int i = 0; i < vertical_profile_points3.size(); i++) {
          out_vel_ver3 << vertical_profile_points3[i][1] << " "
                       << avg_vertical_velocity3[i][0] << " "
                       << avg_vertical_velocity3[i][1] << std::endl;
        }

        output_avg_pressure.close();
        output_avg_stress.close();
        output_Cf.close();
        out_vel_hor.close();
        out_vel_ver1.close();
        out_vel_ver2.close();
        out_vel_ver3.close();
        output_Cp.close();
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }

  }

  template<int dim>
  void NavierStokesProjection<dim>::initialize_points_around_obstacle(const unsigned int n_points, Point<dim> start, double dx){

    obstacle_points.clear();

    double space = 2.0 * dx / (n_points - 1);
    Point<dim> p;

    for(unsigned int i = 0; i < 2*n_points-1; ++i){
      if(i*space < dx){
        p = (dim == 2) ? Point<dim>(start[0], start[1] + i * space) : Point<dim>(start[0], start[1] + i * space, start[2]);
      }
      else if(i*space < 2.0 * dx){
        p = (dim == 2) ? Point<dim>(start[0] + i * space - dx, start[1] + dx) : Point<dim>(start[0] + i * space - dx, start[1] + dx, start[2]);
      }
      else if(i*space < 3.0 * dx){
        p = (dim == 2) ? Point<dim>(start[0] + dx, start[1] + 3.0 * dx - i * space) : Point<dim>(start[0] + dx, start[1] + 3.0 * dx - i * space, start[2]);
      }
      else{
        p = (dim == 2) ? Point<dim>(start[0] + 4.0 * dx - i * space, start[1]) : Point<dim>(start[0] + 4.0 * dx - i * space, start[1], start[2]);
      }
      if(GridTools::find_active_cell_around_point(triangulation, p) != triangulation.end() &&
         GridTools::find_active_cell_around_point(triangulation, p)->is_locally_owned()) {
        obstacle_points.push_back(p);
      }
    }
  }

  template<int dim>
  std::vector<Point<dim>> NavierStokesProjection<dim>::
  initialize_profile_points(double angle, double spacing, Point<dim> start_point, Point<dim> end_point){
    std::vector<Point<dim>> profile_points;
    Point<dim> p = start_point;

    while(p[0] <= end_point[0] && p[1] <= end_point[1]) {
      if(GridTools::find_active_cell_around_point(triangulation, p) != triangulation.end() &&
          GridTools::find_active_cell_around_point(triangulation, p)->is_locally_owned()) {
        profile_points.push_back(p);
      }
      p[0] = p[0] + spacing * std::cos(angle);
      p[1] = p[1] + spacing * std::sin(angle);
    }

    return profile_points;
  }

  template<int dim>
  void NavierStokesProjection<dim>::initialize_nearest_boundary_point_mapping(){
    /*--- Assemble marked vertices of serial triangulation ---*/
    std::vector<bool> global_boundary_vertices(serial_triangulation.n_vertices(), false);
    for (const auto &face : serial_triangulation.active_face_iterators()){
      if(face->at_boundary() && ((no_slip && face->boundary_id() == 3) || face->boundary_id() == 2)){
        for(unsigned int i = 0; i < face->n_vertices(); i++){
          if(!global_boundary_vertices[face->vertex_index(i)] && GridTools::find_active_cell_around_point(triangulation, face->vertex(i)) != triangulation.end() &&
              GridTools::find_active_cell_around_point(triangulation, face->vertex(i))->is_locally_owned()){
            owned_boundary_points[face->vertex_index(i)] = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
          }
          global_boundary_vertices[face->vertex_index(i)] = true;
        }
      }
    }

    /*--- share which boundary points belong to which process ---*/
    auto boundary_points_to_rank_vec = Utilities::MPI::all_gather(MPI_COMM_WORLD, owned_boundary_points);
    for(auto& map : boundary_points_to_rank_vec){
      boundary_points_to_rank.insert(map.begin(), map.end());
    }

    /*--- map each cell to nearest boundary point ---*/
    std::vector<int> boundary_point_indices;
    for(const auto& cell : dof_handler_deltas.active_cell_iterators()) {
      if(cell->is_locally_owned()) {
        int idx = GridTools::find_closest_vertex(serial_triangulation, cell->center(), global_boundary_vertices);

        cell_to_nearest_boundary_point[cell] = idx;
        boundary_point_indices.push_back(idx);
      }
    }

    /*--- Share the map to all other nodes ---*/
    cell_to_nearest_boundary_point_all = Utilities::MPI::all_gather(MPI_COMM_WORLD, boundary_point_indices);
  }

  // pressure average over time
  template<int dim>
  void NavierStokesProjection<dim>::compute_pressure_avg_over_boundary(int n, double height, int n_points) {
    double avg_pres = 0.0;
    for(unsigned int i = 0; i < obstacle_points.size(); i++) {
      for(int j = 0; j < n_points; j++) {
        if(dim == 3)
          obstacle_points[i][2] = j * (height / (n_points-1));
        if(GridTools::find_active_cell_around_point(triangulation, obstacle_points[i]) != triangulation.end() &&
          GridTools::find_active_cell_around_point(triangulation, obstacle_points[i])->is_locally_owned())  {
          avg_pres += VectorTools::point_value(dof_handler_pressure, pres_n, obstacle_points[i]) / n_points;
        }
      }
      if(GridTools::find_active_cell_around_point(triangulation, obstacle_points[i]) != triangulation.end() &&
          GridTools::find_active_cell_around_point(triangulation, obstacle_points[i])->is_locally_owned())  {
        if(n > 1) {
          avg_pressure[i] = ((n-1) * avg_pressure[i] + avg_pres) / n;
        }
        else {
          avg_pressure.push_back(avg_pres);
        }
      }
    }
  }

  // stress average over time over boundary
  template<int dim>
  void NavierStokesProjection<dim>::compute_stress_avg_over_boundary(int n, Point<dim> center, double object_length, double lower_boundary, double upper_boundary) {

    Tensor<1, dim, double> normal_vector;

    for(unsigned int i = 0; i < obstacle_points.size(); i++) {
      if(GridTools::find_active_cell_around_point(triangulation, obstacle_points[i]) != triangulation.end() &&
         GridTools::find_active_cell_around_point(triangulation, obstacle_points[i])->is_locally_owned()) {
        std::vector<Tensor< 1, dim, double >> vel_grad(dim);
        VectorTools::point_gradient(dof_handler_velocity, u_n, obstacle_points[i], vel_grad);

        if(std::abs(obstacle_points[i][1] - lower_boundary) < 1e-10) // south wall
          normal_vector = Tensor<1, dim, double>({0.0, 1.0});
        else if(std::abs(obstacle_points[i][1] - upper_boundary) < 1e-10) // north wall
          normal_vector = Tensor<1, dim, double>({0.0, -1.0});
        else { // square obstacle
          if(obstacle_points[i][0] == center[0] - 0.5*object_length)
            if(obstacle_points[i][1] == center[1] - 0.5*object_length)
              normal_vector = Tensor<1, dim, double>({-0.5*std::sqrt(2.), -0.5*std::sqrt(2.)});
            else if(obstacle_points[i][1] == center[1] + 0.5*object_length)
              normal_vector = Tensor<1, dim, double>({-0.5*std::sqrt(2.), 0.5*std::sqrt(2.)});
            else 
              normal_vector = Tensor<1, dim, double>({-1, 0});
          else if(obstacle_points[i][1] == center[1] - 0.5*object_length)
            if(obstacle_points[i][0] == center[0] + 0.5*object_length)
              normal_vector = Tensor<1, dim, double>({0.5*std::sqrt(2.), -0.5*std::sqrt(2.)});
            else
              normal_vector = Tensor<1, dim, double>({0, -1});
          else if(obstacle_points[i][0] == center[0] + 0.5*object_length)
            if(obstacle_points[i][1] == center[1] + 0.5*object_length)
              normal_vector = Tensor<1, dim, double>({0.5*std::sqrt(2.), 0.5*std::sqrt(2.)});
            else
              normal_vector = Tensor<1, dim, double>({1, 0});
          else if(obstacle_points[i][1] == center[1] + 0.5*object_length)
            normal_vector = Tensor<1, dim, double>({0, 1});
          else
            std::cout << "Error in compute boundary distance for point: " << obstacle_points[i] << std::endl;  
        }
        Tensor< 1, dim, double > tangential_vector = Tensor< 1, dim, double >({normal_vector[1], - normal_vector[0]});

        Tensor< 2, dim, double > vel_grad_tens;
        for(unsigned int i = 0; i < dim; ++i) {
          for(unsigned int j = 0; j < dim; ++j) {
            vel_grad_tens[i][j] = vel_grad[i][j];
          }
        }

        if(n > 1) {
          avg_stress[i] = ((n-1) * avg_stress[i] + vel_grad_tens * normal_vector * tangential_vector) / n;
        }
        else {
          avg_stress.push_back(vel_grad_tens * normal_vector * tangential_vector);
        }
      }
    }
  }

  // velocity average over time
  template<int dim>
  void NavierStokesProjection<dim>::compute_velocity_avg(int n, std::vector<Point<dim>>& points, std::vector<Vector<double>>& avg_velocity) {
    for(unsigned int i = 0; i < points.size(); i++) {
      if(GridTools::find_active_cell_around_point(triangulation, points[i]) != triangulation.end() &&
         GridTools::find_active_cell_around_point(triangulation, points[i])->is_locally_owned()) {

        Vector<double> vel(dim);
        VectorTools::point_value(dof_handler_velocity, u_n, points[i], vel);
        if(n > 1) {
          for(unsigned int d = 0; d < dim; ++d) {
            avg_velocity[i][d] = ((n-1) * avg_velocity[i][d] + vel[d]) / n;
          }
        }
        else {
          avg_velocity.push_back(vel);
        }
      }
    }
  }

  // compute maximal local voriticity
  template<int dim>
  void  NavierStokesProjection<dim>::compute_lipschitz_number() {
    FEValues<dim> fe_values(fe_velocity, quadrature_velocity, update_gradients);
    std::vector<std::vector<Tensor<1, dim, double>>> solution_gradients_velocity(quadrature_velocity.size(),
                                                                                 std::vector<Tensor<1, dim, double>>(dim));

    double max_local_vorticity = std::numeric_limits<double>::min();

    for(const auto& cell: dof_handler_velocity.active_cell_iterators()) {
      if(cell->is_locally_owned()) {
        fe_values.reinit(cell);
        fe_values.get_function_gradients(u_n, solution_gradients_velocity);

        for(unsigned int q = 0; q < quadrature_velocity.size(); ++q)
          max_local_vorticity = std::max(max_local_vorticity,
        std::abs(solution_gradients_velocity[q][1][0] - solution_gradients_velocity[q][0][1])*dt);
      }
    }

    const double lipschitz = Utilities::MPI::max(max_local_vorticity, MPI_COMM_WORLD);
    if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
      output_lipschitz << lipschitz << std::endl;
    }
  }

  // @sect{<code>NavierStokesProjection::compute_lift_and_drag</code>}

  // This routine computes the lift and the drag forces in a non-dimensional framework
  // (so basically for the classical coefficients, it is necessary to multiply by a factor 2).
  //
  template<int dim>
  void NavierStokesProjection<dim>::compute_lift_and_drag() {
    QGauss<dim - 1> face_quadrature_formula(EquationData::degree_p + 2);
    const int n_q_points = face_quadrature_formula.size();

    std::vector<double>                      pressure_values(n_q_points);
    std::vector<std::vector<Tensor<1, dim>>> velocity_gradients(n_q_points, std::vector<Tensor<1, dim>>(dim));

    Tensor<1, dim> normal_vector;
    Tensor<2, dim> fluid_stress;
    Tensor<2, dim> fluid_pressure;
    Tensor<1, dim> forces;

    /*--- We need to compute the integral over the cylinder boundary, so we need to use 'FEFaceValues' instances.
          For the velocity we need the gradients, for the pressure the values. ---*/
    FEFaceValues<dim> fe_face_values_velocity(fe_velocity, face_quadrature_formula,
                                              update_quadrature_points | update_gradients |
                                              update_JxW_values | update_normal_vectors);
    FEFaceValues<dim> fe_face_values_pressure(fe_pressure, face_quadrature_formula, update_values);

    double local_drag = 0.0;
    double local_lift = 0.0;

    /*--- We need to perform a unique loop because the whole stress tensor takes into account contributions of
          velocity and pressure obviously. However, the two dof_handlers are different, so we neede to create an ad-hoc
          iterator for the pressure that we update manually. It is guaranteed that the cells are visited in the same order
          (see the documentation) ---*/
    auto tmp_cell = dof_handler_pressure.begin_active();
    for(const auto& cell : dof_handler_velocity.active_cell_iterators()) {
      if(cell->is_locally_owned()) {
        for(unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face) {
          if(cell->face(face)->at_boundary() && cell->face(face)->boundary_id() == 2) {
            fe_face_values_velocity.reinit(cell, face);
            fe_face_values_pressure.reinit(tmp_cell, face);

            fe_face_values_velocity.get_function_gradients(u_n, velocity_gradients); /*--- velocity gradients ---*/
            fe_face_values_pressure.get_function_values(pres_n, pressure_values); /*--- pressure values ---*/

            for(int q = 0; q < n_q_points; q++) {
              normal_vector = -fe_face_values_velocity.normal_vector(q);

              for(unsigned int d = 0; d < dim; ++ d) {
                fluid_pressure[d][d] = pressure_values[q];
                for(unsigned int k = 0; k < dim; ++k)
                  fluid_stress[d][k] = 1.0/Re*velocity_gradients[q][d][k];
              }
              fluid_stress = fluid_stress - fluid_pressure;

              forces = fluid_stress*normal_vector*fe_face_values_velocity.JxW(q);

              local_drag += forces[0];
              local_lift += forces[1];
            }
          }
        }
      }
      ++tmp_cell;
    }

    /*--- At the end, each processor has computed the contribution to the boundary cells it owns and, therefore,
          we need to sum up all the contributions. ---*/
    const double lift = Utilities::MPI::sum(local_lift, MPI_COMM_WORLD);
    const double drag = Utilities::MPI::sum(local_drag, MPI_COMM_WORLD);
    if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
      output_lift << lift << std::endl;
      output_drag << drag << std::endl;
    }
  }

  template <int dim>
  void NavierStokesProjection<dim>::compute_artificial_force(LinearAlgebra::distributed::Vector<double> & vel) {
    // calculate uniform artificial force
    double force = 0.1 * (1.0 - VectorTools::compute_mean_value(MappingQ1<dim>(), dof_handler_velocity, quadrature_velocity, vel, 0));

    // interpolate force over domain
    VectorTools::interpolate(dof_handler_velocity, Functions::ConstantFunction<dim>(Vector<double>({force, 0.0})), artificial_force);
  }


  template <int dim>
  void NavierStokesProjection<dim>::compute_y_plus(LinearAlgebra::distributed::Vector<double> & vel, double lower_boundary, double upper_boundary, Point<dim> center, double object_length) {
    /*--- calculate gradient of u of all owned boundary points ---*/
    std::map<int, double> boundary_point_to_tau_w;
    for(auto point_pair : owned_boundary_points){
      auto point = serial_triangulation.get_vertices()[point_pair.first];
      std::vector<Tensor< 1, dim, double>> grad_u(dim);
      VectorTools::point_gradient(dof_handler_velocity, vel, point, grad_u);

      // calculate normal vector of boundary point
      Tensor<1, dim, double> normal_vector;
      if(std::abs(point[1] - lower_boundary) < 1e-10) // south wall
        normal_vector = Tensor<1, dim, double>({0.0, 1.0});
      else if(std::abs(point[1] - upper_boundary) < 1e-10) // north wall
        normal_vector = Tensor<1, dim, double>({0.0, -1.0});
      else { // square obstacle
        if(point[0] == center[0] - 0.5*object_length)
          if(point[1] == center[1] - 0.5*object_length)
            normal_vector = Tensor<1, dim, double>({-0.5*std::sqrt(2.), -0.5*std::sqrt(2.)});
          else if(point[1] == center[1] + 0.5*object_length)
            normal_vector = Tensor<1, dim, double>({-0.5*std::sqrt(2.), 0.5*std::sqrt(2.)});
          else 
            normal_vector = Tensor<1, dim, double>({-1, 0});
        else if(point[1] == center[1] - 0.5*object_length)
          if(point[0] == center[0] + 0.5*object_length)
            normal_vector = Tensor<1, dim, double>({0.5*std::sqrt(2.), -0.5*std::sqrt(2.)});
          else
            normal_vector = Tensor<1, dim, double>({0, -1});
        else if(point[0] == center[0] + 0.5*object_length)
          if(point[1] == center[1] + 0.5*object_length)
            normal_vector = Tensor<1, dim, double>({0.5*std::sqrt(2.), 0.5*std::sqrt(2.)});
          else
            normal_vector = Tensor<1, dim, double>({1, 0});
        else if(point[1] == center[1] + 0.5*object_length)
          normal_vector = Tensor<1, dim, double>({0, 1});
        else
          std::cout << "Error in compute boundary distance for point: " << point << std::endl;  
      }

      // calculate tangential vector
      auto tangential_vector = Tensor<1, dim, double>({normal_vector[1], -normal_vector[0]});

      Tensor<1, dim, double> normal_grad_u;
      for(unsigned int idx = 0; idx < grad_u.size(); idx++){
        normal_grad_u[idx] = scalar_product(grad_u[idx], normal_vector);
      }

      boundary_point_to_tau_w[point_pair.first] = std::abs(scalar_product(normal_grad_u, tangential_vector));
    }

    /*--- Share gradient of u to all processes ---*/
    auto boundary_point_to_tau_w_all = Utilities::MPI::all_gather(MPI_COMM_WORLD, boundary_point_to_tau_w);

    /*--- calculate y plus for each owned cell ---*/
    for(const auto& cell : dof_handler_deltas.active_cell_iterators()) {
      if(cell->is_locally_owned()) {
        // calculate distance
        auto boundary_point = serial_triangulation.get_vertices()[cell_to_nearest_boundary_point[cell]];
        double distance_y = boundary_point.distance(cell->center());

        // get grad_u in the normal direction of the boundary point
        double tau_w = boundary_point_to_tau_w_all[boundary_points_to_rank[cell_to_nearest_boundary_point[cell]]][cell_to_nearest_boundary_point[cell]];

        // get dof indices
        std::vector<types::global_dof_index> dof_indices(fe_deltas.dofs_per_cell);
        cell->get_dof_indices(dof_indices);

        // assign y plus
        for(unsigned int idx = 0; idx < dof_indices.size(); ++idx) {
          y_plus(dof_indices[idx]) = std::min(250., distance_y * std::sqrt(Re * tau_w));
          boundary_distance(dof_indices[idx]) = distance_y;
        }
      }
    }
    navier_stokes_matrix.set_y_plus(y_plus);
  }

  // @sect{ <code>NavierStokesProjection::refine_mesh</code>}

  // After finding a good initial guess on the coarse mesh, we hope to
  // decrease the error through refining the mesh. We also need to transfer the current solution to the
  // next mesh using the SolutionTransfer class.
  //
  template <int dim>
  void NavierStokesProjection<dim>::refine_mesh() {
    TimerOutput::Scope t(time_table, "Refine mesh");

    /*--- We first create a proper vector for computing estimator ---*/
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler_velocity, locally_relevant_dofs);
    LinearAlgebra::distributed::Vector<double> tmp_velocity;
    tmp_velocity.reinit(dof_handler_velocity.locally_owned_dofs(), locally_relevant_dofs, MPI_COMM_WORLD);
    tmp_velocity = u_n;
    tmp_velocity.update_ghost_values();

    using Iterator = typename DoFHandler<dim>::active_cell_iterator;
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    /*--- This is basically the indicator per cell computation (see step-50). Since it is not so complciated
          we implement it through a lambda expression ---*/
    const auto cell_worker = [&](const Iterator&   cell,
                                 ScratchData<dim>& scratch_data,
                                 CopyData&         copy_data) {
      FEValues<dim>& fe_values = scratch_data.fe_values; /*--- Here we finally use the 'FEValues' inside ScratchData ---*/
      fe_values.reinit(cell);

      /*--- Compute the gradients for all quadrature points ---*/
      std::vector<std::vector<Tensor<1, dim>>> gradients(fe_values.n_quadrature_points, std::vector<Tensor<1, dim>>(dim));
      fe_values.get_function_gradients(tmp_velocity, gradients);
      copy_data.cell_index = cell->active_cell_index();
      double vorticity_norm_square = 0.0;
      /*--- Loop over quadrature points and evaluate the integral multiplying the vorticty
            by the weights and the determinant of the Jacobian (which are included in 'JxW') ---*/
      for(unsigned k = 0; k < fe_values.n_quadrature_points; ++k) {
        const double vorticity = gradients[k][1][0] - gradients[k][0][1];
        vorticity_norm_square += vorticity*vorticity*fe_values.JxW(k);
      }
      copy_data.value = cell->diameter()*cell->diameter()*vorticity_norm_square;
    };

    const UpdateFlags cell_flags = update_gradients | update_quadrature_points | update_JxW_values;

    const auto copier = [&](const CopyData &copy_data) {
      if(copy_data.cell_index != numbers::invalid_unsigned_int)
        estimated_error_per_cell[copy_data.cell_index] += copy_data.value;
    };

    /*--- Now everything is 'automagically' handled by 'mesh_loop' ---*/
    ScratchData<dim> scratch_data(fe_velocity, EquationData::degree_p + 2, cell_flags);
    CopyData copy_data;
    MeshWorker::mesh_loop(dof_handler_velocity.begin_active(),
                          dof_handler_velocity.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells);

    /*--- Refine grid. In case the refinement level is above a certain value (or the coarsening level is below)
          we clear the flags. ---*/
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(triangulation, estimated_error_per_cell, 0.01, 0.3);
    for(const auto& cell: triangulation.active_cell_iterators()) {
      if(cell->refine_flag_set() && static_cast<unsigned int>(cell->level()) == max_loc_refinements)
        cell->clear_refine_flag();
      if(cell->coarsen_flag_set() && static_cast<unsigned int>(cell->level()) == min_loc_refinements)
        cell->clear_coarsen_flag();
    }
    triangulation.prepare_coarsening_and_refinement();

    /*--- Now we prepare the object for transfering, basically saving the old quantities using SolutionTransfer.
          Since the 'prepare_for_coarsening_and_refinement' method can be called only once, but we have two vectors
          for dof_handler_velocity, we need to put them in an auxiliary vector. ---*/
    std::vector<const LinearAlgebra::distributed::Vector<double>*> velocities;
    velocities.push_back(&u_n);
    velocities.push_back(&u_n_minus_1);
    parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
    solution_transfer_velocity(dof_handler_velocity);
    solution_transfer_velocity.prepare_for_coarsening_and_refinement(velocities);
    parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
    solution_transfer_pressure(dof_handler_pressure);
    solution_transfer_pressure.prepare_for_coarsening_and_refinement(pres_n);

    triangulation.execute_coarsening_and_refinement(); /*--- Effectively perform the remeshing ---*/

    /*--- First DoFHandler objects are set up within the new grid ----*/
    setup_dofs();

    /*--- Interpolate current solutions to new mesh. This is done using auxliary vectors just for safety,
          but the new u_n or pres_n could be used. Again, the only point is that the function 'interpolate'
          can be called once and so the vectors related to 'dof_handler_velocity' have to collected in an auxiliary vector. ---*/
    LinearAlgebra::distributed::Vector<double> transfer_velocity,
                                               transfer_velocity_minus_1,
                                               transfer_pressure;
    transfer_velocity.reinit(u_n);
    transfer_velocity.zero_out_ghost_values();
    transfer_velocity_minus_1.reinit(u_n_minus_1);
    transfer_velocity_minus_1.zero_out_ghost_values();
    transfer_pressure.reinit(pres_n);
    transfer_pressure.zero_out_ghost_values();

    std::vector<LinearAlgebra::distributed::Vector<double>*> transfer_velocities;
    transfer_velocities.push_back(&transfer_velocity);
    transfer_velocities.push_back(&transfer_velocity_minus_1);
    solution_transfer_velocity.interpolate(transfer_velocities);
    transfer_velocity.update_ghost_values();
    transfer_velocity_minus_1.update_ghost_values();
    solution_transfer_pressure.interpolate(transfer_pressure);
    transfer_pressure.update_ghost_values();

    u_n         = transfer_velocity;
    u_n_minus_1 = transfer_velocity_minus_1;
    pres_n      = transfer_pressure;
  }


  // Interpolate the locally refined solution to a mesh with maximal resolution
  // and transfer velocity and pressure.
  //
  template<int dim>
  void NavierStokesProjection<dim>::interpolate_max_res(const unsigned int level) {
    parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
    solution_transfer_velocity(dof_handler_velocity);
    std::vector<const LinearAlgebra::distributed::Vector<double>*> velocities;
    velocities.push_back(&u_n);
    velocities.push_back(&u_n_minus_1);
    solution_transfer_velocity.prepare_for_coarsening_and_refinement(velocities);

    parallel::distributed::SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
    solution_transfer_pressure(dof_handler_pressure);
    solution_transfer_pressure.prepare_for_coarsening_and_refinement(pres_n);

    for(const auto& cell: triangulation.active_cell_iterators_on_level(level)) {
      if(cell->is_locally_owned())
        cell->set_refine_flag();
    }
    triangulation.execute_coarsening_and_refinement();

    setup_dofs();

    LinearAlgebra::distributed::Vector<double> transfer_velocity, transfer_velocity_minus_1,
                                               transfer_pressure;

    transfer_velocity.reinit(u_n);
    transfer_velocity.zero_out_ghost_values();
    transfer_velocity_minus_1.reinit(u_n_minus_1);
    transfer_velocity_minus_1.zero_out_ghost_values();

    transfer_pressure.reinit(pres_n);
    transfer_pressure.zero_out_ghost_values();

    std::vector<LinearAlgebra::distributed::Vector<double>*> transfer_velocities;

    transfer_velocities.push_back(&transfer_velocity);
    transfer_velocities.push_back(&transfer_velocity_minus_1);
    solution_transfer_velocity.interpolate(transfer_velocities);
    transfer_velocity.update_ghost_values();
    transfer_velocity_minus_1.update_ghost_values();

    solution_transfer_pressure.interpolate(transfer_pressure);
    transfer_pressure.update_ghost_values();

    u_n         = transfer_velocity;
    u_n_minus_1 = transfer_velocity_minus_1;
    pres_n      = transfer_pressure;
  }


  // Save maximum resolution to a mesh adapted.
  //
  template<int dim>
  void NavierStokesProjection<dim>::save_max_res() {
    parallel::distributed::Triangulation<dim> triangulation_tmp(MPI_COMM_WORLD);

    if(dim == 2){
      parallel::distributed::Triangulation<dim> tria1(MPI_COMM_WORLD),
                                                  tria2(MPI_COMM_WORLD),
                                                  tria3(MPI_COMM_WORLD),
                                                  tria4(MPI_COMM_WORLD);

      GridGenerator::subdivided_hyper_rectangle(tria1, {60, 19},
                                                Point<dim>(0.0, 0.0),
                                                Point<dim>(30.0, 9.5));
      GridGenerator::subdivided_hyper_rectangle(tria2, {19, 2},
                                                Point<dim>(0.0, 9.5),
                                                Point<dim>(9.5, 10.5));
      GridGenerator::subdivided_hyper_rectangle(tria3, {39, 2},
                                                Point<dim>(10.5, 9.5),
                                                Point<dim>(30.0, 10.5));
      GridGenerator::subdivided_hyper_rectangle(tria4, {60, 19},
                                                Point<dim>(0.0, 10.5),
                                                Point<dim>(30.0, 20.0));
      GridGenerator::merge_triangulations({&tria1, &tria2, &tria3, &tria4},
                                          triangulation, 1e-8, true);

      /*--- Set boundary id ---*/
      for(const auto& face : triangulation.active_face_iterators()) {
        if(face->at_boundary()) {
          const Point<dim> center = face->center();

          // left side
          if(std::abs(center[0] - 0.0) < 1e-10)
            face->set_boundary_id(0);
          // right side
          else if(std::abs(center[0] - 30.0) < 1e-10)
            face->set_boundary_id(1);
          // cylinder boundary
          else if(center[0] < 10.5 + 1e-10 && center[0] > 9.5 - 1e-10 &&
                    center[1] < 10.5 + 1e-10 && center[1] > 9.5 - 1e-10)
            face->set_boundary_id(2);
          // sides of channel
          else {
            Assert(std::abs(center[1] - 0.00) < 1.0e-10 ||
                  std::abs(center[1] - 20.0) < 1.0e-10,
                  ExcInternalError());
            face->set_boundary_id(3);
          }
        }
      }
    }
    else{
      std::cout << "save_max_res() for 3D not implemented" << std::endl;
    }

    triangulation_tmp.refine_global(triangulation.n_global_levels() - 1);

    DoFHandler<dim> dof_handler_velocity_tmp(triangulation_tmp);
    DoFHandler<dim> dof_handler_pressure_tmp(triangulation_tmp);
    dof_handler_velocity_tmp.distribute_dofs(fe_velocity);
    dof_handler_pressure_tmp.distribute_dofs(fe_pressure);

    LinearAlgebra::distributed::Vector<double> u_n_tmp,
                                               pres_n_tmp;
    u_n_tmp.reinit(dof_handler_velocity_tmp.n_dofs());
    pres_n_tmp.reinit(dof_handler_pressure_tmp.n_dofs());

    DataOut<dim> data_out;
    std::vector<std::string> velocity_names(dim, "v");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation_velocity(dim, DataComponentInterpretation::component_is_part_of_vector);
    VectorTools::interpolate_to_different_mesh(dof_handler_velocity, u_n, dof_handler_velocity_tmp, u_n_tmp);
    u_n_tmp.update_ghost_values();
    data_out.add_data_vector(dof_handler_velocity_tmp, u_n_tmp, velocity_names, component_interpretation_velocity);
    VectorTools::interpolate_to_different_mesh(dof_handler_pressure, pres_n, dof_handler_pressure_tmp, pres_n_tmp);
    pres_n_tmp.update_ghost_values();
    data_out.add_data_vector(dof_handler_pressure_tmp, pres_n_tmp, "p", {DataComponentInterpretation::component_is_scalar});
    PostprocessorVorticity<dim> postprocessor;
    data_out.add_data_vector(dof_handler_velocity_tmp, u_n_tmp, postprocessor);

    data_out.build_patches(MappingQ<dim>(EquationData::degree_p + 1, false),
			   EquationData::degree_p + 1, DataOut<dim>::curved_inner_cells);
    const std::string output = "./" + saving_dir + "/solution_max_res_end.vtu";
    data_out.write_vtu_in_parallel(output, MPI_COMM_WORLD);
  }

  template<int dim>
  void read_statistics_velocity(std::vector<Point<dim>>& points, int axis, std::vector<Vector<double>>& values, std::string filename) {
    std::ifstream infile;
    double idx;
    Vector<double> val(dim);

    values.resize(points.size());

    infile.open(filename);
    while (infile >> idx){
      infile >> val[0] >> val[1];
      unsigned int it = std::find_if(points.begin(), points.end(),
        [&](Point<dim> p){ return std::abs(p[axis] - idx) < 1.0e-10; }) - points.begin();

      if(it != points.size()){
        values[it] = val;
      }
    }
    infile.close();
  }

  template<int dim>
  void read_statistics(std::vector<Point<dim>> & points, std::vector<double> & values, std::string filename) {
    std::ifstream infile;
    Point<dim> p;
    double val;

    values.resize(points.size());

    infile.open(filename);
    while (infile >> p[0] && infile >> p[1]){
      infile >> val;
      unsigned int it = std::find_if(points.begin(), points.end(),
        [&](Point<dim> point){ return p == point; }) - points.begin();

      if(it != points.size()){
        values[it] = val;
      }
      else
        std::cout << "Error in read_statistics" << std::endl;
    }
    infile.close();
  }

  // @sect{ <code>NavierStokesProjection::run</code> }

  // This is the time marching function, which starting at <code>t_0</code>
  // advances in time using the projection method with time step <code>dt</code>
  // until <code>T</code>.
  //
  // Its second parameter, <code>verbose</code> indicates whether the function
  // should output information what it is doing at any given moment:
  // we use the ConditionalOStream class to do that for us.
  //
  template<int dim>
  void NavierStokesProjection<dim>::run(const bool verbose, const unsigned int output_interval) {
    ConditionalOStream verbose_cout(std::cout, verbose && Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
    
    initialize_nearest_boundary_point_mapping();

    Point<dim> center;
    double radius, length, height, y_start;

    radius = 0.5;
    height = 20.0;
    length = 30.0; 

    if(import_mesh) {
      center =  Point<dim>(0.0, 0.0);
      y_start = -10.0;
    }
    else {
      center =  Point<dim>(10.0, 10.0);
      y_start = 0.0;
    }

    // verbose_cout << " initialize statistics points" << std::endl;
    // initialize_points_around_obstacle(200, Point<dim>(center[0] - radius, center[1] - radius), 2.0 * radius);
    // horizontal_wake_points   = initialize_profile_points(0.0, 0.01, Point<dim>(center[0] + radius, 0.5 * height), Point<dim>(length, 0.5 * height));
    // vertical_profile_points1 = initialize_profile_points(0.5 * numbers::PI, 0.01, Point<dim>(center[0] + 1.05 * 2.0 * radius, 0.0), Point<dim>(center[1] + 1.05 * 2.0 * radius, height));
    // vertical_profile_points2 = initialize_profile_points(0.5 * numbers::PI, 0.01, Point<dim>(center[0] + 1.54 * 2.0 * radius, 0.0), Point<dim>(center[1] + 1.54 * 2.0 * radius, height));
    // vertical_profile_points3 = initialize_profile_points(0.5 * numbers::PI, 0.01, Point<dim>(center[0] + 2.02 * 2.0 * radius, 0.0), Point<dim>(center[1] + 2.02 * 2.0 * radius, height));

    double time = t_0 + dt;
    unsigned int n = 1;
    if(restart && !as_initial_conditions) {
      n    = step_restart;
      time = time_restart;

      // read_statistics_velocity(horizontal_wake_points, 0, avg_horizontal_velocity, "./" + saving_dir + "/out_vel_hor.dat");
      // read_statistics_velocity(vertical_profile_points1, 1, avg_vertical_velocity1, "./" + saving_dir + "/out_vel_ver1.dat");
      // read_statistics_velocity(vertical_profile_points2, 1, avg_vertical_velocity2, "./" + saving_dir + "/out_vel_ver2.dat");
      // read_statistics_velocity(vertical_profile_points3, 1, avg_vertical_velocity3, "./" + saving_dir + "/out_vel_ver3.dat");

      // read_statistics(obstacle_points, avg_stress, "./" + saving_dir + "/avg_stress.dat");
      // read_statistics(obstacle_points, avg_pressure, "./" + saving_dir + "/avg_p.dat");
    }
    else {
      compute_y_plus(u_n, y_start, y_start + height, center, 2.0 * radius);

      output_results(1);

      // compute_pressure_avg_over_boundary(n);
      // compute_stress_avg_over_boundary(n, center, 2.0 * radius);
      // compute_velocity_avg(n, horizontal_wake_points, avg_horizontal_velocity);
      // compute_velocity_avg(n, vertical_profile_points1, avg_vertical_velocity1);
      // compute_velocity_avg(n, vertical_profile_points2, avg_vertical_velocity2);
      // compute_velocity_avg(n, vertical_profile_points3, avg_vertical_velocity3);
    }
    while(std::abs(T - time) > 1e-10) {
      time += dt;
      n++;
      pcout << "Step = " << n << " Time = " << time << std::endl;

      /*--- compute artificial force for first step ---*/
      compute_artificial_force(u_n);

      /*--- First stage of TR-BDF2 and we start by setting the proper flag ---*/
      TR_BDF2_stage = 1;
      navier_stokes_matrix.set_TR_BDF2_stage(TR_BDF2_stage);
      for(unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
        mg_matrices[level].set_TR_BDF2_stage(TR_BDF2_stage);

      verbose_cout << "  Interpolating the velocity stage 1" << std::endl;
      // auto start = std::chrono::high_resolution_clock::now();
      interpolate_velocity();
      // auto stop = std::chrono::high_resolution_clock::now();
      // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
      // compute y+
      verbose_cout << "  Computing y+ stage 1" << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      compute_y_plus(u_n, y_start, y_start + height, center, 2.0 * radius);
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;

      verbose_cout << "  Diffusion Step stage 1 " << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      diffusion_step();
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      verbose_cout << "  Projection Step stage 1" << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      project_grad(1);
      u_tmp.equ(gamma*dt, u_tmp);
      u_star += u_tmp; /*--- In the rhs of the projection step we need u_star + gamma*dt*grad(pres_n) and we save it into u_star ---*/
      projection_step();
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      verbose_cout << "  Updating the Velocity stage 1" << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      u_n_gamma.equ(1.0, u_star);
      project_grad(2);
      grad_pres_int.equ(1.0, u_tmp); /*--- We save grad(pres_int), because we will need it soon ---*/
      u_tmp.equ(-gamma*dt, u_tmp);
      u_n_gamma += u_tmp; /*--- u_n_gamma = u_star - gamma*dt*grad(pres_int) ---*/
      u_n_minus_1 = u_n;
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      /*--- Second stage of TR-BDF2 ---*/
      TR_BDF2_stage = 2;
      for(unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
        mg_matrices[level].set_TR_BDF2_stage(TR_BDF2_stage);
      navier_stokes_matrix.set_TR_BDF2_stage(TR_BDF2_stage);

      verbose_cout << "  Interpolating the velocity stage 2" << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      interpolate_velocity();
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      // compute y+
      verbose_cout << "  Computing y+ stage 2" << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      compute_y_plus(u_n_gamma, y_start, y_start + height, center, 2.0 * radius);
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      verbose_cout << "  Diffusion Step stage 2 " << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      diffusion_step();
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      verbose_cout << "  Projection Step stage 2" << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      u_tmp.equ((1.0 - gamma)*dt, grad_pres_int);
      u_star += u_tmp;  /*--- In the rhs of the projection step we need u_star + (1 - gamma)*dt*grad(pres_int) ---*/
      projection_step();
      // stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      verbose_cout << "  Updating the Velocity stage 2" << std::endl;
      // start = std::chrono::high_resolution_clock::now();
      u_n.equ(1.0, u_star);
      project_grad(1);
      u_tmp.equ((gamma - 1.0)*dt, u_tmp);
      u_n += u_tmp;  /*--- u_n = u_star - (1 - gamma)*dt*grad(pres_n) ---*/
      // stop = std::chrono::high_resolution_clock::now();

      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
      // pcout << duration << std::endl;
 
      const double max_vel = get_maximal_velocity();
      pcout<< "Maximal velocity = " << max_vel << std::endl;
      /*--- The Courant number is computed taking into account the polynomial degree for the velocity ---*/
      pcout << "CFL = " << dt*max_vel*(EquationData::degree_p + 1)*
                           std::sqrt(dim)/GridTools::minimal_cell_diameter(triangulation) << std::endl;

      compute_lift_and_drag();

      // // compute time average of parameters along different points
      // compute_pressure_avg_over_boundary(n);
      // compute_stress_avg_over_boundary(n, center, 2.0 * radius);
      // compute_velocity_avg(n, horizontal_wake_points, avg_horizontal_velocity);
      // compute_velocity_avg(n, vertical_profile_points1, avg_vertical_velocity1);
      // compute_velocity_avg(n, vertical_profile_points2, avg_vertical_velocity2);
      // compute_velocity_avg(n, vertical_profile_points3, avg_vertical_velocity3);

      // // compute lipschitz number at every timestep
      // compute_lipschitz_number();

      if(n % output_interval == 0) {
        verbose_cout << "Plotting Solution final" << std::endl;
        output_results(n);
        // output_statistics(center);
      }
      /*--- In case dt is not a multiple of T, we reduce dt in order to end up at T ---*/
      if(T - time < dt && T - time > 1e-10) {
        dt = T - time;
        navier_stokes_matrix.set_dt(dt);
        for(unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
          mg_matrices[level].set_dt(dt);
      }
      /*--- Perform the refinement if desired ---*/
      if(refinement_iterations > 0 && n % refinement_iterations == 0) {
        verbose_cout << "Refining mesh" << std::endl;
        refine_mesh();
      }
      /*--- Modify the Reynolds number if desired ---*/
      if(modify_Reynolds) {
        if(restart && time <= time_restart + 0.05 + 1e-10) {
          const double Re_tmp = (time - time_restart)/(0.05)*(3800.0) + 100.0;
          navier_stokes_matrix.set_Reynolds(Re_tmp);
          verbose_cout << "Reynolds set to " << Re_tmp << std::endl;
        }
        if(!restart && time <= t_0 + dt + 0.05 + 1e-10) {
          const double Re_tmp = (time - t_0 - dt)/(0.05)*(3800.0) + 100.0;
          navier_stokes_matrix.set_Reynolds(Re_tmp);
          verbose_cout << "Reynolds set to " << Re_tmp << std::endl;
        }
      }
    }

    if(n % output_interval != 0) {
      verbose_cout << "Plotting Solution final" << std::endl;
      output_results(n);
    }
    if(refinement_iterations > 0) {
      for(unsigned int lev = 0; lev < triangulation.n_global_levels() - 1; ++ lev)
        interpolate_max_res(lev);
      save_max_res();
    }
  }

} // namespace NS_TRBDF2


// @sect{ The main function }

// The main function looks very much like in all the other tutorial programs. We first initialize MPI,
// we initialize the class 'NavierStokesProjection' with the dimension as template parameter and then
// let the method 'run' do the job.
//
int main(int argc, char *argv[]) {
  try {
    using namespace NS_TRBDF2;

    RunTimeParameters::Data_Storage data;
    data.read_data("parameter-file.prm");

    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, -1);

    const auto& curr_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    deallog.depth_console(data.verbose && curr_rank == 0 ? 2 : 0);

    NavierStokesProjection<2> test(data);
    test.run(data.verbose, data.output_interval);

    if(curr_rank == 0)
      std::cout << "----------------------------------------------------"
                << std::endl
                << "Apparently everything went fine!" << std::endl
                << "Don't forget to brush your teeth :-)" << std::endl
                << std::endl;

    return 0;
  }
  catch(std::exception &exc) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

}
