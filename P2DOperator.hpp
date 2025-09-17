#include "mfem.hpp"
#include "equations/ElectrolytePotential.hpp"
#include "equations/ElectrolyteConcentration.hpp"
#include "equations/SolidPotential.hpp"
#include "equations/SolidConcentration.hpp"
#include "coefficients/ExchangeCurrentCoefficient.hpp"
#include "coefficients/OpenCircuitPotentialCoefficient.hpp"
#include "coefficients/OverPotentialCoefficient.hpp"
#include "coefficients/ReactionCurrentCoefficient.hpp"
#include "RegionalCurrent.hpp"

using namespace mfem;

#pragma once

class P2DOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace * &_x_fespace;
   Array<ParFiniteElementSpace *> &_r_fespace;

   ElectrolytePotential     *  _ep;
   SolidPotential           *  _sp;
   ElectrolyteConcentration *  _ec;
   Array<SolidConcentration *> _sc;

   /// Gridfunctions defined over x_fespace (3 macro eqs plus _surface_ concentration)
   ParGridFunction * _ep_gf;
   ParGridFunction * _sp_gf;
   ParGridFunction * _ec_gf;
   ParGridFunction * _sc_gf;

   /// Big-enough array for the surface concentrations of the two SPM(e) particles
   Array<real_t> _sc_array{/* UNKNOWN */ 0., /* NE */ 0., /* SEP */ 0., /* PE */ 0.};

   /// Coefficients for derived, i.e. not solved for, quantities
   ReactionCurrentCoefficient * _j;
   ExchangeCurrentCoefficient * _jex;
   OpenCircuitPotentialCoefficient * _ocp;
   OverPotentialCoefficient * _op;

   /// For the four gridfunctions over _x_fespace (3 macro eqs plus _surface_ concentration)
   Array<int> _block_offsets;
   /// For solution true vector (3 macros eqs plus NPAR _radial_ concentrations)
   Array<int> _block_trueOffsets;
   /// For rhs true vectors (2 macro eqs)
   Array<int> _potential_trueOffsets;
   /// For rhs true vectors (1 macro eq plus NPAR _radial_ concentrations)
   Array<int> _concentration_trueOffsets;

   /// System matrices for concentration and potential eqs
   BlockOperator * _Ac, * _Ap;

   /// Block vector for the dofs of quantities defined over _x_fespace (3 macro eqs plus _surface_ concentration)
   BlockVector _l;

   /// Reference to solution true dof vector
   BlockVector & _x;

   /// Reference to current time
   real_t & _t;

   /// Reference to current timestep
   real_t & _dt;

   /// Reference to the ODE solver we use to integrate in time
   ODESolver & _ode_solver;

   /// Implicit solver for T = M + dt K
   CGSolver _Solver;
   /// Preconditioner for the implicit solver
   HypreSmoother _Prec;

   /// Auxiliary rhs vectors for concentrations and potential eqs
   mutable BlockVector _bc, _bp;

   /// File to write temporary data to
   std::ofstream _file;

public:
   P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace, const unsigned &ndofs,
               BlockVector &x, real_t & t, real_t & dt, ODESolver & ode_solver);

   virtual void Mult(const Vector &x, Vector &dx_dt) const override {};

   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;

   void Step();
   void SetGridFunctionsFromTrueVectors();
   void SetSurfaceConcentration();

   virtual void UpdatePotentialEquations();
   virtual void UpdateConcentrationEquations();

   real_t GetElectrodeReactionCurrent(const Region &r, const int &sign);
   Array<real_t> GetParticleReactionCurrent();

   /// Construct coefficients for derived quantities
   void ConstructReactionCurrent();
   void ConstructExchangeCurrent();
   void ConstructOpenCircuitPotential();
   void ConstructOverPotential();

   /// SPM(e) helpers
   const real_t & GetSurfaceConcentration(const Region &r);
   const real_t & GetReactionCurrent(const Region &r);
   const real_t & GetExchangeCurrent(const Region &r);
   const real_t & GetOpenCircuitPotential(const Region &r);
   const real_t & GetOverPotential(const Region &r);

   real_t GetVoltage();

   virtual void GetParticleDofs(Array<int> & particle_dofs, Array<Region> & particle_regions, Array<int> & particle_offsets);

   virtual ~P2DOperator()
   {
      delete _ep;
      delete _sp;
      delete _ec;
      for (unsigned p = 0; p < NPAR; p++)
         delete _sc[p];
      delete _ep_gf;
      delete _sp_gf;
      delete _ec_gf;
      delete _sc_gf;
      delete _Ac;
      delete _Ap;
      _file.close();
   }
};
