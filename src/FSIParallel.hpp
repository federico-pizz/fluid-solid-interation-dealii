// Base
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

// Distributed
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/grid_refinement.h>

// DoFs
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

// Finite Elements
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

// Preconditioning
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_ilu.h>

// Grid
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

// hp
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

// Linear Algebra
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/vector.h>

// Numerics
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

// Trilinos
#include <deal.II/base/mpi.h>
#include <deal.II/lac/trilinos_precondition.h>

// Standard library
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sched.h>
#include <thread>
#include <unistd.h>

// Timing
#include <deal.II/base/timer.h>

using namespace dealii;

class FSI {
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

  class InletVelocity : public Function<dim> {
  public:
    InletVelocity() : Function<dim>(dim + 1 + dim) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component) const override {
      Assert(component < this->n_components,
             ExcIndexRange(component, 0, this->n_components));

      if (component == dim - 1)
        switch (dim) {
        case 2:
          return std::sin(numbers::PI * p[0]);
        case 3:
          return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]);
        default:
          Assert(false, ExcNotImplemented());
        }
      return 0;
    }

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override {
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = this->value(p, c);
    }
  };

  class FSIPreconditioner {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void initialize(
        const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
        const TrilinosWrappers::SparseMatrix &pressure_mass_,
        const TrilinosWrappers::SparseMatrix &displacement_stiffness_,
        const TrilinosWrappers::SparseMatrix &B10_,
        const TrilinosWrappers::SparseMatrix &B20_,
        const TrilinosWrappers::SparseMatrix &B21_,
        const std::vector<std::vector<bool>> &displacement_constant_modes) {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass = &pressure_mass_;
      B10 = &B10_;
      B20 = &B20_;
      B21 = &B21_;
      displacement_stiffness = &displacement_stiffness_;

      
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data_pres;
      amg_data_pres.elliptic = true;
      amg_data_pres.higher_order_elements = false;
      amg_data_pres.aggregation_threshold = 1e-2;
      preconditioner_pressure.initialize(pressure_mass_, amg_data_pres);

      /**
       * The AMG takes some options for optimization.
       * When setting aggregation_threshold (to 0.02 for example)
       * we specify that all entries that are more than two percent
       * of size of some diagonal pivots in that row should form one
       * coarse grid point. This parameter can be fine-tuned to gain
       * performance. A rule of thumb is that larger values will decrease
       * the number of iterations, but increase the costs per iteration.
       */
      // Displacement AMG
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
      amg_data.constant_modes = displacement_constant_modes;
      amg_data.elliptic = true;
      amg_data.higher_order_elements = true;
      amg_data.smoother_sweeps = 3;
      amg_data.w_cycle = false;
      amg_data.aggregation_threshold = 1e-3;
      preconditioner_displacement.initialize(displacement_stiffness_, amg_data);

      //preconditioner_velocity.initialize(velocity_stiffness_);
      
      // Velocity AMG
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data_vel;
      amg_data_vel.elliptic = true;
      amg_data_vel.higher_order_elements = true;
      amg_data_vel.smoother_sweeps = 2;
      amg_data_vel.aggregation_threshold = 1e-2;
      preconditioner_velocity.initialize(velocity_stiffness_, amg_data_vel);
      
    }

    // Application of the preconditioner.
    void vmult(TrilinosWrappers::MPI::BlockVector &dst,
               const TrilinosWrappers::MPI::BlockVector &src) const {
      tmp_p.reinit(src.block(1));
      tmp_d.reinit(src.block(2));
      intermediate_tmp.reinit(src.block(2));
      {
        // Solve fluid system
        SolverControl solver_control_vel(2000, 1e-2 * src.block(0).l2_norm());

        SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_vel(
            solver_control_vel);
        dst.block(0) = 0;

        solver_cg_vel.solve(*velocity_stiffness, dst.block(0), src.block(0),
                               preconditioner_velocity);
      }

      // Compute Pressure residual: r_p = src[1] - B10 * u
      B10->vmult(tmp_p, dst.block(0));
      tmp_p.sadd(-1.0, src.block(1));

      {
        SolverControl solver_control_pres(2000, 1e-2 * tmp_p.l2_norm());

        SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pres(
            solver_control_pres);
        dst.block(1) = 0;

        solver_cg_pres.solve(*pressure_mass, dst.block(1), tmp_p,
                             preconditioner_pressure);
      }

      B20->vmult(tmp_d, dst.block(0));            // fluid velocity contribution
      B21->vmult(intermediate_tmp, dst.block(1)); // fluid pressure contribution

      tmp_d.add(1.0, intermediate_tmp); // Sum total coupling effect (B * px)
      tmp_d.sadd(-1.0, src.block(2));   // residual = y - B*px

      {
        SolverControl solver_control_disp(2000, 1e-2 * tmp_d.l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_disp(
            solver_control_disp);
        dst.block(2) = 0;

        solver_gmres_disp.solve(*displacement_stiffness, dst.block(2), tmp_d,
                                preconditioner_displacement);
      }

      // preconditioner_displacement.vmult(dst.block(2), tmp_d);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Displacement stiffness matrix
    const TrilinosWrappers::SparseMatrix *displacement_stiffness;
    const TrilinosWrappers::SparseMatrix *B10;
    const TrilinosWrappers::SparseMatrix *B20;
    const TrilinosWrappers::SparseMatrix *B21;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionAMG preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionAMG preconditioner_pressure;

    // Preconditioner used for the solid block
    TrilinosWrappers::PreconditionAMG preconditioner_displacement;

    // Temporary vector.
    mutable TrilinosWrappers::MPI::Vector tmp_p;
    mutable TrilinosWrappers::MPI::Vector tmp_d;
    mutable TrilinosWrappers::MPI::Vector intermediate_tmp;
  };

  // Constructor
  FSI(const unsigned int &degree_velocity_,
      const unsigned int &degree_pressure_,
      const unsigned int &degree_displacement_, const double &nu_,
      const double &mu_, const double &lambda_)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0), nu(nu_), mu(mu_), lambda(lambda_),
        degree_velocity(degree_velocity_), degree_pressure(degree_pressure_),
        degree_displacement(degree_displacement_), mesh(MPI_COMM_WORLD) {}

  // Setup system (mesh, FE space, DoF handler, and linear system).
  void setup();

  // Assemble system matrix and right-hand side vector.
  // We also assemble the pressure mass matrix (needed for the preconditioner).
  void assemble_system();

  // Assemble the interface terms between fluid and solid domains.
  void assemble_interface_term(
      const FEFaceValuesBase<dim> &elasticity_fe_face_values,
      const FEFaceValuesBase<dim> &stokes_fe_face_values,
      std::vector<Tensor<1, dim>> &elasticity_phi,
      std::vector<Tensor<2, dim>> &stokes_grad_phi_u,
      std::vector<double> &stokes_phi_p,
      FullMatrix<double> &local_interface_matrix) const;

  // Solve the linear system.
  void solve();

  // Write solution to output file.
  void output(const unsigned int refinement_cycle);

  // The function that generates the mesh
  void make_grid();

  // Refine mesh based on error estimation.
  void refine_mesh();

  // Main solver loop.
  void run();

protected:
  // MPI /////////////////////////////////////////////////////////////////////
  const unsigned int mpi_size;
  const unsigned int mpi_rank;
  ConditionalOStream pcout;

  // Domain identifiers for material ID tagging.
  static constexpr types::material_id fluid_domain_id = 0;
  static constexpr types::material_id solid_domain_id = 1;

  // Physical and material parameters. ////////////////////////////////////////

  // Kinematic viscosity [m2/s]. - deal-ii uses 2 for some reason
  const double nu;

  

  // Lamé parameters.
  const double mu;
  const double lambda;

  // Forcing term.
  Tensor<1, dim> f;

  // Dirichlet datum. Can be omitted in our problem
  // FunctionG function_g;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Discretization parameters and objects. ///////////////////////////////////

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Polynomial degree used for displacement
  const unsigned int degree_displacement;

  // Mesh.
  parallel::distributed::Triangulation<dim> mesh;

  // Finite element spaces.
  std::unique_ptr<FiniteElement<dim>> stokes_fe;
  std::unique_ptr<FiniteElement<dim>> elasticity_fe;

  hp::FECollection<dim> fe_collection;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> stokes_quadrature;
  std::unique_ptr<Quadrature<dim>> elasticity_quadrature;

  // Quadrature formula for face integrals.
  std::unique_ptr<Quadrature<dim - 1>> common_face_quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // Sparsity pattern
  BlockSparsityPattern sparsity_pattern;
  BlockSparsityPattern pressure_mass_sparsity;

  // Constraints object for handling boundary conditions and hanging nodes
  // (check what these are).
  AffineConstraints<double> constraints;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix for the Stokes preconditioner
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;
  TrilinosWrappers::MPI::BlockVector solution_owned;
  TrilinosWrappers::MPI::BlockVector solution;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  std::vector<IndexSet> block_owned_dofs;
  std::vector<IndexSet> block_relevant_dofs;

  // Helper functions.
  // ///////////////////////////////////////////////////////////
  static bool
  cell_is_in_fluid_domain(const typename DoFHandler<dim>::cell_iterator &cell);

  static bool
  cell_is_in_solid_domain(const typename DoFHandler<dim>::cell_iterator &cell);
};
