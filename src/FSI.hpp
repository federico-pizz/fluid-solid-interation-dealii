// Base
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

// Distributed
#include <deal.II/distributed/fully_distributed_tria.h>

// DoFs
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

// Finite Elements
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

// Grid
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

// hp
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

// Linear Algebra
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
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

// Standard library
#include <fstream>
#include <iostream>

using namespace dealii;

class FSI
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;  
  
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
      : Function<dim>(dim + 1 + dim)
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component) const override
    {
      Assert(component < this->n_components, 
          ExcIndexRange(component, 0, this->n_components));
      
      if (component == dim - 1)
        switch (dim)
          {
            case 2:
              return std::sin(numbers::PI * p[0]);
            case 3:
              return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]);
            default:
              Assert(false, ExcNotImplemented());
      }
      return 0;
    }

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = this->value(p, c);
    }
  };


  // Constructor
  // do we want to use template like deal-ii?
  FSI(const std::string &mesh_file_name_,
      const unsigned int &degree_velocity_,
      const unsigned int &degree_pressure_,
      const unsigned int &degree_displacement_)
  : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
  , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
  , pcout(std::cout, mpi_rank == 0)
  , mesh_file_name(mesh_file_name_)
  , degree_velocity(degree_velocity_)
  , degree_pressure(degree_pressure_)
  , degree_displacement(degree_displacement_)
  , mesh(MPI_COMM_WORLD)
  {}  

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

  // Refine mesh based on error estimation.
  void refine_mesh();

  // Main solver loop.
  void run();

protected:
  // Domain identifiers for material ID tagging.
  static constexpr types::material_id fluid_domain_id = 0;
  static constexpr types::material_id solid_domain_id = 1;

  // MPI information. ////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Physical and material parameters. ////////////////////////////////////////

  // Kinematic viscosity [m2/s]. - deal-ii uses 2 for some reason
  const double nu = 2;

  // Outlet pressure [Pa].
  const double p_out = 10;

  // Lamé parameters.
  const double mu     = 1.0;
  const double lambda = 10.0;

  // Forcing term.
  Tensor<1, dim> f;

  // Dirichlet datum. Can be omitted in our problem
  // FunctionG function_g; 

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Discretization parameters and objects. ///////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Polynomial degree used for displacement
  const unsigned int degree_displacement;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

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

  // Constraints object for handling boundary conditions and hanging nodes (check what these are).
  AffineConstraints<double> constraints;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // System matrix (FSI operator)
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Mass matrix (velocity + displacement only)
  TrilinosWrappers::BlockSparseMatrix mass_matrix;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  // Helper functions. ///////////////////////////////////////////////////////////

  static bool cell_is_in_fluid_domain(
    const typename DoFHandler<dim>::cell_iterator &cell);

  static bool cell_is_in_solid_domain(
    const typename DoFHandler<dim>::cell_iterator &cell);

private:
  void assemble_mass_matrix();
};


