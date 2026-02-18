#include "FSI.hpp"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // Problem parameters
  // Physical compile-time dimension and mesh size
  constexpr unsigned int dim = 2;
  const unsigned int n_el = 8;

  const unsigned int degree_velocity     = 2;
  const unsigned int degree_pressure     = 1;
  const unsigned int degree_displacement = 2;

  // Physical parameters
  const double nu     = 2.0;  // Kinematic viscosity
  const double p_out  = 10.0; // Outlet pressure
  const double mu     = 1.0;  // Lamé mu
  const double lambda = 10.0; // Lamé lambda

  // Create FSI problem
  FSI<dim> problem(n_el,
                   degree_velocity,
                   degree_pressure,
                   degree_displacement,
                   nu,
                   p_out,
                   mu,
                   lambda);

  problem.run();

  return 0;
}