#include "P2DOperator.hpp"

P2DOperator::P2DOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace, const unsigned &ndofs,
                         BlockVector &x, real_t & t, real_t & dt, ODESolver & ode_solver)
   : TimeDependentOperator(ndofs, (real_t) 0.0), _x_fespace(x_fespace), _r_fespace(r_fespace), _Ac(NULL), _Ap(NULL),
     _x(x), _t(t), _dt(dt), _ode_solver(ode_solver), _Solver(_x_fespace->GetComm()), _file("data.csv")
{
   const real_t rel_tol = 1e-16;

   _Solver.iterative_mode = false;
   _Solver.SetRelTol(rel_tol);
   _Solver.SetAbsTol(0.0);
   _Solver.SetMaxIter(100);
   _Solver.SetPrintLevel(0);
   //_Solver.SetPreconditioner(_Prec);

   _block_offsets.SetSize(NMACRO + 1 + 1);
   _block_trueOffsets.SetSize(NEQS + 1);
   _potential_trueOffsets.SetSize(NMACROP + 1);
   _concentration_trueOffsets.SetSize(NMACROC + NPAR + 1);

   _block_offsets[0] = 0;
   _block_offsets[EP + 1] = _x_fespace->GetVSize();
   _block_offsets[SP + 1] = _x_fespace->GetVSize();
   _block_offsets[EC + 1] = _x_fespace->GetVSize();
   _block_offsets[SC + 1] = _x_fespace->GetVSize();

   _block_trueOffsets[0] = 0;
   _block_trueOffsets[EP + 1] = _x_fespace->GetTrueVSize();
   _block_trueOffsets[SP + 1] = _x_fespace->GetTrueVSize();
   _block_trueOffsets[EC + 1] = _x_fespace->GetTrueVSize();

   _potential_trueOffsets[0] = 0;
   _potential_trueOffsets[EPP + 1] = _x_fespace->GetTrueVSize();
   _potential_trueOffsets[SPP + 1] = _x_fespace->GetTrueVSize();

   _concentration_trueOffsets[0] = 0;
   _concentration_trueOffsets[ECC + 1] = _x_fespace->GetTrueVSize();

   for (unsigned p = 0; p < NPAR; p++)
   {
      _block_trueOffsets[SC + 1 + p] = _r_fespace[p]->GetTrueVSize();
      _concentration_trueOffsets[SCC + 1 + p] = _r_fespace[p]->GetTrueVSize();
   }

   _block_offsets.PartialSum();
   _block_trueOffsets.PartialSum();
   _potential_trueOffsets.PartialSum();
   _concentration_trueOffsets.PartialSum();

   if (!Mpi::WorldRank())
   {
      std::cout << "Variables: " << NEQS << std::endl;
      std::cout << "Unknowns (rank 0): " << _block_trueOffsets[NEQS] << std::endl;
   }

   // Set offsets for full dof vector
   _l.Update(_block_offsets); _l = 0.;

   // Initialise gridfunctions to use the appropriate section of the full dof vector _l
   _ep_gf = new ParGridFunction(_x_fespace, _l, _block_offsets[EP]);
   _sp_gf = new ParGridFunction(_x_fespace, _l, _block_offsets[SP]);
   _ec_gf = new ParGridFunction(_x_fespace, _l, _block_offsets[EC]);
   _sc_gf = new ParGridFunction(_x_fespace, _l, _block_offsets[SC]);

   // Set offsets for solution and rhs (potential and concentration) true vectors
   _x.Update(_block_trueOffsets); _x = 0.;
   _bp.Update(_potential_trueOffsets);
   _bc.Update(_concentration_trueOffsets);

   // Construct equation ojects, first the 3 macro equations, then the NPAR micro eqs
   _ep = new ElectrolytePotential(*_x_fespace);
   _sp = new SolidPotential(*_x_fespace);
   _ec = new ElectrolyteConcentration(*_x_fespace);

   if (M == SPM || M == SPMe)
   {
      _sc.Append(new SolidConcentration(*_r_fespace[0], 0, 0, -1, NE));
      _sc.Append(new SolidConcentration(*_r_fespace[1], 1, 0, -1, PE));
   }
   else
   {
      Array<int> particle_dofs, particle_offsets;
      Array<Region> particle_regions;
      GetParticleDofs(particle_dofs, particle_regions, particle_offsets);

      for (unsigned p = 0; p < NPAR; p++)
      {
         auto rank_iter = std::upper_bound(particle_offsets.begin(), particle_offsets.end(), p);
         unsigned rank = std::distance(particle_offsets.begin(), rank_iter) - 1;
         bool owned = rank == Mpi::WorldRank();

         unsigned offset = particle_offsets[Mpi::WorldRank()];
         int dof = owned ? particle_dofs[p - offset] : -1;
         Region region = owned ? particle_regions[p - offset] : UNKNOWN;

         _sc.Append(new SolidConcentration(*_r_fespace[p], p, rank, dof, region));
      }
   }

   ConstructExchangeCurrent();
   ConstructOpenCircuitPotential();
   ConstructOverPotential();
   ConstructReactionCurrent();
}

void P2DOperator::ImplicitSolve(const real_t dt,
                                const Vector &x, Vector &dx_dt)
{
   // Solve the equation:
   //   M dx_dt = -K(x + dt*dx_dt) <=> (M + dt K) dx_dt = -Kx
   // for dx_dt, where K is linearized by using x from the previous timestep

   // Logically split dx_dt in two parts: potentials and concentrations
   Array<int> offsets({_block_trueOffsets[P], _block_trueOffsets[C], _block_trueOffsets[NEQS]});
   BlockVector dx_dt_blocked(dx_dt, offsets);
   Vector & dxp_dt(dx_dt_blocked.GetBlock(0));
   Vector & dxc_dt(dx_dt_blocked.GetBlock(1));

   if (M == P2D)
   {
      dx_dt = 0.;

      // for (SCL)
      {
         // assemble each individual block of _Ap and _bp
         _x.Add(_dt, dx_dt);
         SetGridFunctionsFromTrueVectors();
         // compute absolute potentials in dedicated member function
         UpdatePotentialEquations();

         real_t Inp = GetElectrodeReactionCurrent(NE,  1.0);
         real_t Inn = GetElectrodeReactionCurrent(NE,  1.0);
         real_t Ipp = GetElectrodeReactionCurrent(PE, -1.0);
         real_t Ipn = GetElectrodeReactionCurrent(PE, -1.0);

         real_t phi_eo = -2.0 * log(I + sqrt(4.0 * Inp * Inn + I * I)/(2.0 * Inp));
         real_t phi_po = phi_eo + 2.0 * log(-I + sqrt(4.0 * Ipp * Ipn + I * I)/(2.0 * Ipp));

         _x.Add(-_dt, dx_dt);

         // put _Ap and _bp together
         _Ap->SetDiagonalBlock(EPP, new HypreParMatrix(_ep->GetK()));
         _Ap->SetDiagonalBlock(SPP, new HypreParMatrix(_sp->GetK()));
         _bp.GetBlock(EPP) = _ep->GetZ();
         _bp.GetBlock(SPP) = _sp->GetZ();

         // solve for dxp_dt (potentials rate)
         _Solver.SetOperator(*_Ap);
         _Solver.Mult(_bp, dxp_dt);
      }
   }

   // assemble each individual block of _Ac and _bc
   UpdateConcentrationEquations();

   // put _Ac and _bc together
   _Ac->SetDiagonalBlock(ECC, Add(1, _ec->GetM(), dt, _ec->GetK()));
   _bc.GetBlock(ECC) = _ec->GetZ();
   for (unsigned p = 0; p < NPAR; p++)
   {
      _Ac->SetDiagonalBlock(SCC + p, Add(1, _sc[p]->GetM(), dt, _sc[p]->GetK()));
      _bc.GetBlock(SCC + p) = _sc[p]->GetZ();
   }

   // solve for dxc_dt (concentrations rate)
   _Solver.SetOperator(*_Ac);
   _Solver.Mult(_bc, dxc_dt);
}

void P2DOperator::Step()
{
   _ode_solver.Step(_x, _t, _dt);
   SetGridFunctionsFromTrueVectors();
}

void P2DOperator::SetGridFunctionsFromTrueVectors()
{
   _ep_gf->SetFromTrueDofs(_x.GetBlock(EP));
   _sp_gf->SetFromTrueDofs(_x.GetBlock(SP));
   _ec_gf->SetFromTrueDofs(_x.GetBlock(EC));
   SetSurfaceConcentration();
}

void P2DOperator::UpdatePotentialEquations()
{
   // Rebuild _Ap, destroys owned (i.e. all) blocks
   delete _Ap;
   _Ap = new BlockOperator(_potential_trueOffsets);
   _Ap->owns_blocks = 1;

   _ep->Update(_x, *_j, _dt);
   _sp->Update(_x, *_j, _dt);
}

void P2DOperator::UpdateConcentrationEquations()
{
   // Rebuild _Ac, destroys owned (i.e. all) blocks
   delete _Ac;
   _Ac = new BlockOperator(_concentration_trueOffsets);
   _Ac->owns_blocks = 1;

   _ec->Update(_x, *_j, _dt);
   const Array<real_t> & j = GetParticleReactionCurrent();
   for (unsigned p = 0; p < NPAR; p++)
      _sc[p]->Update(_x, ConstantCoefficient(j[p]), _dt);
}

//
// Surface Concentration
//

const real_t & P2DOperator::GetSurfaceConcentration(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get constant surface concentration, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot get constant surface concentration, only negative (NE) and positive electrodes (PE) are supported.");
   return _sc_array[r];
}

void P2DOperator::SetSurfaceConcentration()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         for (unsigned p = 0; p < NPAR; p++)
         {
            Region r = _sc[p]->GetParticleRegion();
            real_t sc = (r == NE ? CN0 : CP0) + _sc[p]->SurfaceConcentration(_x);
            _sc_array[r] = sc;
         }
         break;
      case P2D:
         for (unsigned p = 0; p < NPAR; p++)
         {
            Region r = _sc[p]->GetParticleRegion();
            real_t sc = (r == NE ? CN0 : CP0) + _sc[p]->SurfaceConcentration(_x);
            if (_sc[p]->IsParticleOwned())
               (*_sc_gf)(_sc[p]->GetParticleDof()) = sc;
         }
         // Apply prolongation after restriction. Might be unnecessary, but guarantees
         // all processors have the right information for all their local dofs.
         _sc_gf->SetFromTrueVector();
         break;
   }
}

//
// Partial Reaction Current for each electrode
//

real_t P2DOperator::GetElectrodeReactionCurrent(const Region &r, const int &sign)
{
   MFEM_ASSERT(r == NE || r == PE, "Cannot get partial electrode reaction current, only negative (NE) and positive electrodes (PE) are supported.");

   Vector amask({r == NE ? AN : 0., 0., r == PE ? AP : 0.});
   PWConstCoefficient a(amask);
   ProductCoefficient ajex(a, *_jex);
   TransformedCoefficient I(&ajex, _op, [=](real_t ajex, real_t op) { return ajex * exp( sign * 0.5 * op ); });
   QuadratureSpace x_qspace(_x_fespace->GetParMesh(), 2 * _x_fespace->FEColl()->GetOrder());
   return x_qspace.Integrate(I);
}

//
// Reaction Current for each particle
//

Array<real_t> P2DOperator::GetParticleReactionCurrent()
{
   Array<real_t> j(NPAR);

   switch (M)
   {
      case SPM:
      case SPMe:
         for (unsigned p = 0; p < NPAR; p++)
            j[p] = GetReactionCurrent(_sc[p]->GetParticleRegion());
         break;
      case P2D:
         ParGridFunction j_gf(_x_fespace);
         j_gf.ProjectCoefficient(*_j);

         for (unsigned p = 0; p < NPAR; p++)
         {
            j[p] = _sc[p]->IsParticleOwned() ? j_gf(_sc[p]->GetParticleDof()) : 0;

            if (_sc[p]->GetParticleRank() == _sc[p]->GetSurfaceRank())
               continue;

            MPI_Request request;
            if (_sc[p]->IsParticleOwned())
               MPI_Isend(&j[p], 1, MFEM_MPI_REAL_T, _sc[p]->GetSurfaceRank(), 1, MPI_COMM_WORLD, &request);

            if (_sc[p]->IsSurfaceOwned())
               MPI_Recv(&j[p], 1, MFEM_MPI_REAL_T, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (_sc[p]->IsParticleOwned())
               MPI_Wait(&request, MPI_STATUS_IGNORE);
         }
         break;
   }

   return j;
}

//
// Reaction Current
//

const real_t & P2DOperator::GetReactionCurrent(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant reaction current, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant reaction current, only negative (NE) and positive electrodes (PE) are supported.")
   return _j->Eval()(r);
}

void P2DOperator::ConstructReactionCurrent()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         _j = new ReactionCurrentCoefficient();
         break;
      case P2D:
         _j = new ReactionCurrentCoefficient(T, *_jex, *_op); // TODO: add absolute potentials as inputs
         break;
   }
}

//
// Exchange Current
//

const real_t & P2DOperator::GetExchangeCurrent(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant exchange current, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant exchange current, only negative (NE) and positive electrodes (PE) are supported.");
   return _jex->Eval()(r);
}

void P2DOperator::ConstructExchangeCurrent()
{
   switch (M)
   {
      case SPM:
         _jex = new ExchangeCurrentCoefficient(KN, KP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE), CE0);
         break;
      case SPMe:
         _jex = new ExchangeCurrentCoefficient(KN, KP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE), *_ec_gf);
         break;
      case P2D:
         _jex = new ExchangeCurrentCoefficient(KN, KP, *_sc_gf, *_ec_gf);
         break;
   }
}

//
// Open Circuit Potential
//

const real_t & P2DOperator::GetOpenCircuitPotential(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant open circuit potential, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant open circuit potential, only negative (NE) and positive electrodes (PE) are supported.")
   return _ocp->Eval()(r);
}

void P2DOperator::ConstructOpenCircuitPotential()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         _ocp = new OpenCircuitPotentialCoefficient(UN, UP, GetSurfaceConcentration(NE), GetSurfaceConcentration(PE));
         break;
      case P2D:
         _ocp = new OpenCircuitPotentialCoefficient(UN, UP, *_sc_gf);
         break;
   }
}

//
// OverPotential
//

const real_t & P2DOperator::GetOverPotential(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot calculate constant  overpotential, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot calculate constant overpotential, only negative (NE) and positive electrodes (PE) are supported.")
   return _op->Eval()(r);
}

void P2DOperator::ConstructOverPotential()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         _op = new OverPotentialCoefficient(T, *_jex);
         break;
      case P2D:
         _op = new OverPotentialCoefficient(*_sp_gf, *_ep_gf, *_ocp);
         break;
   }
}

//
// Voltage
//

real_t P2DOperator::GetVoltage()
{
   real_t voltage;

   switch (M)
   {
      case SPM:
      case SPMe:
         // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
         voltage = GetOpenCircuitPotential(PE) - GetOpenCircuitPotential(NE) +
                   (GetOverPotential(PE) - GetOverPotential(NE)) * phi_scale;
         break;
      case P2D:
         voltage = 0.;
         break;
   }

   if (Mpi::Root())
   {
      std::cout << "[Rank " << Mpi::WorldRank() << "]"
                << " Voltage = " << voltage << std::endl;

      // Print file headings first time function is called.
      static bool writeFileHeadings = true;
      if (writeFileHeadings) {
         _file << "t"  << ", \t"
               << "cn" << ", \t"
               << "cp" << ", \t"
               << "voltage"
               << std::endl;

         writeFileHeadings = false;
      }

      // Print data to file.
      _file << _t                          << ", \t"
            << GetSurfaceConcentration(NE) << ", \t"
            << GetSurfaceConcentration(PE) << ", \t"
            << voltage
            << std::endl;
   }

   _sc[1]->DebuggingCheck(_x);

   return voltage;
}

void P2DOperator::GetParticleDofs(Array<int> & particle_dofs, Array<Region> & particle_regions, Array<int> & particle_offsets)
{
   std::set<std::pair<int, Region>> electrode_dofs_set;
   std::set<int> sep_gdofs_set;
   for (int e = 0; e < _x_fespace->GetNE(); e++)
   {
      Array<int> dofs;
      _x_fespace->GetElementDofs(e, dofs);
      switch (Region r = Region(_x_fespace->GetAttribute(e)))
      {
         case NE:
         case PE:
            for (int d: dofs)
               electrode_dofs_set.insert({d, r});
            break;
         case SEP:
            for (int d: dofs)
               sep_gdofs_set.insert(_x_fespace->GetGlobalTDofNumber(d));
            break;
      }
   }

   unsigned max_sep_gdofs = NSEP * _x_fespace->FEColl()->GetOrder() + 1;
   Array<int> sep_gdofs(max_sep_gdofs); sep_gdofs = -1;
   std::copy(sep_gdofs_set.begin(), sep_gdofs_set.end(), sep_gdofs.begin());

   Array<int> all_sep_gdofs(max_sep_gdofs * Mpi::WorldSize());
   MPI_Allgather(sep_gdofs.GetData(), max_sep_gdofs, MPI_INT,
                 all_sep_gdofs.GetData(), max_sep_gdofs, MPI_INT, MPI_COMM_WORLD);

   Array<int> boundary_dofs;
   _x_fespace->GetBoundaryTrueDofs(boundary_dofs);

   for (auto [dof, region]: electrode_dofs_set)
   {
      int ltdof = _x_fespace->GetLocalTDofNumber(dof);
      int gtdof = _x_fespace->GetGlobalTDofNumber(dof);
      if (ltdof != -1 && boundary_dofs.Find(ltdof) == -1 && all_sep_gdofs.Find(gtdof) == -1)
      {
         particle_dofs.Append(dof);
         particle_regions.Append(region);
      }
   }

   int my_particles = particle_dofs.Size();
   particle_offsets.SetSize(Mpi::WorldSize());
   MPI_Allgather(&my_particles, 1, MPI_INT, particle_offsets.GetData(), 1, MPI_INT, MPI_COMM_WORLD);

   particle_offsets.Prepend(0);
   particle_offsets.PartialSum();
   MFEM_ASSERT(particle_offsets[Mpi::WorldSize()] == NPAR, "Failed to distribute particles across processors.");
}
