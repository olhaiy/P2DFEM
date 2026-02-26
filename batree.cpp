// Batree - An MFEM-based SPM, SPMe and P2D solver

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "operators/EChemOperator.hpp"

using namespace mfem;

int
main(int argc, char * argv[])
{
  // Initialize MPI and HYPRE.
  Mpi::Init(argc, argv);
  Hypre::Init();

  // Parse command-line options.
  std::string model = "SPM";
  int order = 1;
  int ode_solver_type = 1;
  real_t t_final = 3600.0;
  real_t dt = 1.0;
  int output_steps = 5;

  int precision = 8;
  std::cout.precision(precision);

  OptionsParser args(argc, argv);
  args.AddOption(&model, "-m", "--model", "Electrochemical model: SPM, SPMe, or P2D.");
  args.AddOption(&order, "-o", "--order", "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type,
                 "-s",
                 "--ode-solver",
                 "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&output_steps, "-os", "--output-steps", "Output every n-th timestep.");
  args.Parse();
  if (!args.Good())
  {
    args.PrintUsage(std::cout);
    return 1;
  }

  if (Mpi::Root())
    args.PrintOptions(std::cout);

  // Define the ODE solver used for time integration. Several implicit
  // singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
  // explicit Runge-Kutta methods are available in MFEM. For now, we only
  // support implicit methods and we have only tested Backward Euler.
  ODESolver * ode_solver;
  switch (ode_solver_type)
  {
    // Implicit L-stable methods
    case 1:
      ode_solver = new BackwardEulerSolver;
      break;
    case 2:
      ode_solver = new SDIRK23Solver(2);
      break;
    case 3:
      ode_solver = new SDIRK33Solver;
      break;
    // Explicit methods
    case 11:
      ode_solver = new ForwardEulerSolver;
      break;
    case 12:
      ode_solver = new RK2Solver(0.5);
      break; // midpoint method
    case 13:
      ode_solver = new RK3SSPSolver;
      break;
    case 14:
      ode_solver = new RK4Solver;
      break;
    case 15:
      ode_solver = new GeneralizedAlphaSolver(0.5);
      break;
    // Implicit A-stable methods (not L-stable)
    case 22:
      ode_solver = new ImplicitMidpointSolver;
      break;
    case 23:
      ode_solver = new SDIRK23Solver;
      break;
    case 24:
      ode_solver = new SDIRK34Solver;
      break;
    default:
      std::cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
      return 1;
  }

  // Initialise grid and layout properties dependent on the electrochemical model and FE order
  init_settings(model, order);

  // Build the 1d mesh for the macro problem and tag its elements according to their region.
  // Define the parallel mesh by a partitioning of the serial mesh.
  // Once the parallel mesh is defined, the serial mesh can be deleted.
  Mesh x_smesh = Mesh::MakeCartesian1D(NX);
  for (unsigned i = 0; i < NX; i++)
    x_smesh.SetAttribute(i, i < NNE ? NE : i < NNE + NSEP ? SEP : PE);
  ParMesh * x_pmesh = new ParMesh(MPI_COMM_WORLD, x_smesh);
  x_smesh.Clear(); // the serial mesh is no longer needed

  // Build one 1d mesh for each particle, i.e. for each of the micro problems.
  // Define each parallel mesh by a partitioning of the respective serial mesh.
  // Once each parallel mesh is defined, the respective serial mesh can be deleted.
  Array<ParMesh *> r_pmesh(NPAR);
  for (unsigned p = 0; p < NPAR; p++)
  {
    Mesh r_smesh = Mesh::MakeCartesian1D(NR);
    r_pmesh[p] = new ParMesh(MPI_COMM_WORLD, r_smesh);
    r_smesh.Clear(); // the serial mesh is no longer needed
  }

  // Define the H1 finite element spaces representing concentrations/potentials
  H1_FECollection fe_coll(order, /*dim*/ 1);
  ParFiniteElementSpace * x_h1space = new ParFiniteElementSpace(x_pmesh, &fe_coll);
  Array<ParFiniteElementSpace *> r_h1space(NPAR);
  for (unsigned p = 0; p < NPAR; p++)
    r_h1space[p] = new ParFiniteElementSpace(r_pmesh[p], &fe_coll);

  // Get the total number of dofs in the system (including boundaries), for
  // both the macro and micro problems, across all processors. This is for
  // reporting purposes only.
  {
    HYPRE_BigInt fe_size_global = NMACRO * x_h1space->GlobalTrueVSize();
    for (unsigned p = 0; p < NPAR; p++)
      fe_size_global += r_h1space[p]->GlobalTrueVSize();

    if (Mpi::Root())
      std::cout << "Unknowns (total): " << fe_size_global << std::endl;
  }

  // Get the number of dofs in the system (including boundaries), for
  // both the macro and micro problems, _owned_ by this processor.
  HYPRE_BigInt fe_size_owned = NMACRO * x_h1space->GetTrueVSize();
  for (unsigned p = 0; p < NPAR; p++)
    fe_size_owned += r_h1space[p]->GetTrueVSize();

  // Initialize the ElectroChemistry operator.
  real_t t = 0.0;
  BlockVector x;
  EChemOperator oper(x_h1space, r_h1space, fe_size_owned, x, t, dt, *ode_solver);

  // Perform time-integration (looping over the time iterations, ti, with a
  // time-step dt).
  ode_solver->Init(oper);

  bool last_step = false;
  for (int ti = 1; !last_step; ti++)
  {
    last_step = t + dt >= t_final - dt / 2;

    oper.Step();
    real_t V = oper.GetVoltage();
    // TODO: Stop sim at cutoff voltage

    if (Mpi::Root() && output_steps && ti == 1)
      std::cout << "step\ttime[s]\tvoltage[V]" << std::endl;

    if (Mpi::Root() && output_steps && (last_step || (ti % output_steps) == 0))
      std::cout << ti << "\t" << t << "\t" << V << std::endl;
  }

  // Free the used memory.
  delete ode_solver;
  delete x_pmesh;
  for (unsigned p = 0; p < NPAR; p++)
    delete r_pmesh[p];

  return 0;
}
