#include "FSIParallel.hpp"

int main(int argc, char *argv[]) {

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  // Problem parameters
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  const unsigned int degree_displacement = 2;

  // Physical parameters
  const double nu = 2.0;      // Kinematic viscosity
  const double mu = 1.0;      // Lamé mu
  const double lambda = 10.0; // Lamé lambda

  // Create FSI problem
  FSI problem(degree_velocity,
              degree_pressure,
              degree_displacement,
              nu,
              mu,
              lambda);

  // Run the simulation
  problem.run();

  return 0;
}