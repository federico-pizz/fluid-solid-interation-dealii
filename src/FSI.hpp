#ifndef FSI_HPP
#define FSI_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

using namespace dealii;

class FSI {
public:

  // Dirichlet Boundary Condition for the Inlet Velocity.
  class InletVelocity : public Function<dim> {
  public:
    InletVelocity() : Function<dim>(dim + 1 + dim) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component) const override {
      if (component == dim - 1) {
        switch (dim) {
        case 2:
          return std::sin(numbers::PI * p[0]);
        case 3:
          return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]);
        default:
          Assert(false, ExcNotImplemented());
        }
      }
      return 0;
    }

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override {
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = this->value(p, c);
    }
  };

  // Preconditioner class. RE-DO
  class FSIPreconditioner {
  public:
    void initialize(
        const TrilinosWrappers::BlockSparseMatrix &system_matrix,
        const TrilinosWrappers::BlockSparseMatrix &pressure_mass,
        const std::vector<std::vector<bool>> &displacement_constant_modes) {
      velocity_stiffness = &system_matrix.block(0, 0);
      this->pressure_mass = &pressure_mass.block(1, 1);
      displacement_stiffness = &system_matrix.block(2, 2);
      B10 = &system_matrix.block(1, 0);
      B20 = &system_matrix.block(2, 0);
      B21 = &system_matrix.block(2, 1);

      // Trilinos AMG preconditioners
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
      amg_data.constant_modes = displacement_constant_modes;
      amg_data.elliptic = true;
      amg_data.higher_order_elements = false;
      amg_data.smoother_sweeps = 3;
      amg_data.w_cycle = false;
      amg_data.aggregation_threshold = 1e-3;
      // Use the displacement block (2,2) for AMG initialization
      preconditioner_displacement.initialize(*displacement_stiffness, amg_data);

      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data_vel;
      amg_data_vel.elliptic = true;
      amg_data_vel.higher_order_elements = true;
      amg_data_vel.smoother_sweeps = 2;
      amg_data_vel.aggregation_threshold = 1e-3;
      // Use the velocity block (0,0) for AMG initialization
      preconditioner_velocity.initialize(*velocity_stiffness, amg_data_vel);

      // Use the pressure mass matrix (1,1) for Jacobi initialization
      preconditioner_pressure.initialize(*this->pressure_mass);
    }

    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const {
      // Velocity block
      SolverControl solver_control_vel(2000, 1e-2 * src.block(0).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_vel(solver_control_vel);
      dst.block(0) = 0;
      solver_cg_vel.solve(*velocity_stiffness, dst.block(0), src.block(0),
                          preconditioner_velocity);

      // Intermediate step
      tmp_p.reinit(src.block(1));
      B10->vmult(tmp_p, dst.block(0));
      tmp_p.sadd(-1.0, src.block(1));

      // Pressure block
      SolverControl solver_control_pres(2000, 1e-2 * tmp_p.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pres(
          solver_control_pres);
      dst.block(1) = 0;
      solver_cg_pres.solve(*pressure_mass, dst.block(1), tmp_p,
                           preconditioner_pressure);

      // Intermediate step
      tmp_d.reinit(src.block(2));
      intermediate_tmp.reinit(src.block(2));
      B20->vmult(tmp_d, dst.block(0));
      B21->vmult(intermediate_tmp, dst.block(1));
      tmp_d.sadd(1.0, intermediate_tmp);
      tmp_d.sadd(-1.0, src.block(2));

      // Displacement block
      SolverControl solver_control_disp(2000, 1e-2 * tmp_d.l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_disp(
          solver_control_disp);
      dst.block(2) = 0;
      solver_gmres_disp.solve(*displacement_stiffness, dst.block(2), tmp_d,
                              preconditioner_displacement);
    }

  protected:
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;
    const TrilinosWrappers::SparseMatrix *displacement_stiffness;
    const TrilinosWrappers::SparseMatrix *B10;
    const TrilinosWrappers::SparseMatrix *B20;
    const TrilinosWrappers::SparseMatrix *B21;
    TrilinosWrappers::PreconditionAMG preconditioner_velocity;
    TrilinosWrappers::PreconditionAMG preconditioner_displacement;
    TrilinosWrappers::PreconditionJacobi preconditioner_pressure;
    const TrilinosWrappers::SparseMatrix *pressure_mass;
    mutable TrilinosWrappers::MPI::Vector tmp_p;
    mutable TrilinosWrappers::MPI::Vector tmp_d;
    mutable TrilinosWrappers::MPI::Vector intermediate_tmp;
  };

  // Constructor with mesh file.
  FSI(const std::string  &mesh_file_name_,
      const unsigned int  n_el_,
      const unsigned int &degree_velocity_,
      const unsigned int &degree_pressure_,
      const unsigned int &degree_displacement_,
      const double       &nu_,
      const double       &p_out_,
      const double       &mu_,
      const double       &lambda_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name(mesh_file_name_)
    , n_el(n_el_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , degree_displacement(degree_displacement_)
    , nu(nu_)
    , p_out(p_out_)
    , mu(mu_)
    , lambda(lambda_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Constructor with generated grid.
  FSI(const unsigned int  n_el_,
      const unsigned int &degree_velocity_,
      const unsigned int &degree_pressure_,
      const unsigned int &degree_displacement_,
      const double       &nu_,
      const double       &p_out_,
      const double       &mu_,
      const double       &lambda_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name("")
    , n_el(n_el_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , degree_displacement(degree_displacement_)
    , nu(nu_)
    , p_out(p_out_)
    , mu(mu_)
    , lambda(lambda_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Main run loop.
  void run();

protected:
  // Initialization.
  void setup();

  void make_grid();

  void refine_mesh();

  void assemble_system();

  void solve();

  void output(const unsigned int refinement_cycle);

  // Determine physics for each cell.
  static bool
  cell_is_in_fluid_domain(const typename DoFHandler<dim>::cell_iterator &cell);

  static bool
  cell_is_in_solid_domain(const typename DoFHandler<dim>::cell_iterator &cell);

  // Domain identifiers.
  static constexpr types::material_id fluid_domain_id = 0;
  static constexpr types::material_id solid_domain_id = 1;

  // MPI information. //////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Physical and material parameters. /////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  const double nu;

  // Outlet pressure [Pa].
  const double p_out;

  // Lamé parameters.
  const double mu;
  const double lambda;

  // Forcing term.
  Tensor<1, dim> f;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Discretization parameters. ////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Number of elements per direction for generated grid.
  const unsigned int n_el;

  // Polynomial degrees.
  const unsigned int degree_velocity;
  const unsigned int degree_pressure;
  const unsigned int degree_displacement;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // DoF handler and FE collections.
  DoFHandler<dim> dof_handler;
  hp::FECollection<dim> fe_collection;
  std::unique_ptr<FESystem<dim>> stokes_fe;
  std::unique_ptr<FESystem<dim>> elasticity_fe;
  std::unique_ptr<QGauss<dim>> stokes_quadrature;
  std::unique_ptr<QGauss<dim>> elasticity_quadrature;
  std::unique_ptr<QGauss<dim - 1>> common_face_quadrature;

  // Constraints and index sets.
  AffineConstraints<double> constraints;
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  std::vector<IndexSet> block_owned_dofs;
  std::vector<IndexSet> block_relevant_dofs;

  // Linear system. ////////////////////////////////////////////////////////////

  // System matrix (FSI operator).
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix (needed for preconditioning).
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Mass matrix (velocity + displacement only).
  TrilinosWrappers::BlockSparseMatrix mass_matrix;

  // Right-hand side vector.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  // Preconditioner.
  FSIPreconditioner preconditioner;

private:
  // Internal helpers.
  void assemble_interface_term(
      const FEFaceValuesBase<dim> &elasticity_fe_face_values,
      const FEFaceValuesBase<dim> &stokes_fe_face_values,
      std::vector<Tensor<1, dim>> &elasticity_phi,
      std::vector<Tensor<2, dim>> &stokes_grad_phi_u,
      std::vector<double> &stokes_phi_p,
      FullMatrix<double> &local_interface_matrix) const;
};

#endif /* FSI_HPP */
