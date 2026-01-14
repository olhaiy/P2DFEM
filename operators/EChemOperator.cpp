#include "operators/EChemOperator.hpp"

EChemOperator::EChemOperator(ParFiniteElementSpace * &x_fespace, Array<ParFiniteElementSpace *> &r_fespace, const unsigned &ndofs,
                             BlockVector &x, real_t & t, real_t & dt, ODESolver & ode_solver)
   : TimeDependentOperator(ndofs, (real_t) 0.0), _x_fespace(x_fespace), _r_fespace(r_fespace), _Ac(NULL), _Ap(NULL),
     _x(x), _t(t), _dt(dt), _ode_solver(ode_solver), _Solver(_x_fespace->GetComm()), _file("data.csv")
{
   const real_t rel_tol = 1e-16;

   _Solver.iterative_mode = false;
   _Solver.SetRelTol(rel_tol);
   _Solver.SetAbsTol(0.0);
   _Solver.SetMaxIter(500);
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
   _x.Update(_block_trueOffsets); _x = 0.; _x.GetBlock(EC) = CE0;
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

void EChemOperator::ImplicitSolve(const real_t dt,
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
   static unsigned iter = 0;

   if (M == P2D)
   {
      dx_dt = 0.;
      ParGridFunction j_gf(_x_fespace);
      ConstantCoefficient zero(0.);

      IntegrationRule ir(4);
      for (size_t i = 0; i < ir.GetNPoints(); i++)
         ir.IntPoint(i).weight = 1.;
      const IntegrationRule * irs[Geometry::Type::NUM_GEOMETRIES];
      irs[Geometry::Type::SEGMENT] = &ir;

      do
      {
         std::cout << "SCL: " << ++iter << std::endl;
         j_gf.ProjectCoefficient(*_jex);
         std::cout << "jex_gf" << std::endl; j_gf.Print();

         // save previous iteration reaction current
         j_gf.ProjectCoefficient(*_j);
         std::cout << "j_gf" << std::endl; j_gf.Print();

         // restore solution true dof vector; this is paramount to guarantee we
         // use the old timestep in the rhs (which appears as a consequence of
         // our rate formulation); note, however, that any coefficients
         // dependent on the potentials should use the gridfunctions set below,
         // which are on the new timestep, just like _j, NOT any gridfunctions
         // obtained afresh from true vector _x which is now on the old timestep
         _x.Add(-_dt, dx_dt);

         // assemble each individual block of _Ap and _bp
         UpdatePotentialEquations();

         // put _Ap and _bp together
         _Ap->SetDiagonalBlock(EPP, new HypreParMatrix(_ep->GetK()));
         _Ap->SetDiagonalBlock(SPP, new HypreParMatrix(_sp->GetK()));
         _bp.GetBlock(EPP) = _ep->GetZ();
         _bp.GetBlock(SPP) = _sp->GetZ();

         // solve for dxp_dt (potentials rate)
         _Solver.SetOperator(*_Ap);
         _Solver.Mult(_bp, dxp_dt);

         // temporarily advance solution true dof vector to set gridfunctions
         _x.Add(_dt, dx_dt);
         std::cout << "SP" << std::endl; _x.GetBlock(SP).Print();
         std::cout << "EP" << std::endl;  _x.GetBlock(EP).Print();
         SetGridFunctionsFromTrueVectors();

         std::cout << "Error: " << std::setprecision(16) << std::fixed << j_gf.ComputeL2Error(*_j, irs) << std::endl;
         std::cout << "Norm: " << std::setprecision(16) << std::fixed << j_gf.ComputeL2Error(zero, irs) << std::endl;
      }
      while (j_gf.ComputeL2Error(*_j, irs) > _threshold * j_gf.ComputeL2Error(zero, irs));

      // restore solution true dof vector
      _x.Add(-_dt, dx_dt);
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

void EChemOperator::Step()
{
   _ode_solver.Step(_x, _t, _dt);
   SetGridFunctionsFromTrueVectors();
}

void EChemOperator::SetGridFunctionsFromTrueVectors()
{
   _ep_gf->SetFromTrueDofs(_x.GetBlock(EP));
   _sp_gf->SetFromTrueDofs(_x.GetBlock(SP));
   _ec_gf->SetFromTrueDofs(_x.GetBlock(EC));
   SetSurfaceConcentration();
   std::cout << "Surf concentration: " << std::endl;
   _sc_gf->Print();
   SetReferencePotential();
}

void EChemOperator::UpdatePotentialEquations()
{
   // Rebuild _Ap, destroys owned (i.e. all) blocks
   delete _Ap;
   _Ap = new BlockOperator(_potential_trueOffsets);
   _Ap->owns_blocks = 1;

   _ep->Update(_x, *_j, _dt);
   _sp->Update(_x, *_j, _dt);
}

void EChemOperator::UpdateConcentrationEquations()
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

const real_t & EChemOperator::GetSurfaceConcentration(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get surface concentration, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot get surface concentration, only negative (NE) and positive electrodes (PE) are supported.");
   return _sc_array[r];
}

void EChemOperator::SetSurfaceConcentration()
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
// Reference Potential
//

const real_t & EChemOperator::GetReferencePotential(const Region &r)
{
   MFEM_ASSERT(M == P2D, "Cannot get reference potential, only P2D is supported.");
   MFEM_ASSERT(r == E || r == NE || r == PE, "Cannot get reference potential, only electrolyte (E), negative (NE) and positive electrodes (PE) are supported.");
   return _rp_array[r];
}

void EChemOperator::SetReferencePotential()
{
    switch (M)
   {
      case SPM:
      case SPMe:
         break;
      case P2D:
         _rp_array[E] = 0.;
         _rp_array[PE] = 0.;

         real_t Inp = GetElectrodeReactionCurrent(NE,  1.0); std::cout << Inp << std::endl;
         real_t Inn = GetElectrodeReactionCurrent(NE, -1.0); std::cout << Inn << std::endl;
         real_t Ipp = GetElectrodeReactionCurrent(PE,  1.0); std::cout << Ipp << std::endl;
         real_t Ipn = GetElectrodeReactionCurrent(PE, -1.0); std::cout << Ipn << std::endl;

         _rp_array[E] = -2.0 * T * log(( I + sqrt(4.0 * Inp * Inn + I * I))/(2.0 * Inp));
         _rp_array[PE] = 2.0 * T * log((-I + sqrt(4.0 * Ipp * Ipn + I * I))/(2.0 * Ipp)) + _rp_array[E];
         std::cout << "Ref pots: " << _rp_array[E] << " " << _rp_array[PE] << std::endl;
         break;
   }
}

//
// Partial Reaction Current for each electrode
//

real_t EChemOperator::GetElectrodeReactionCurrent(const Region &r, const int &sign)
{
   MFEM_ASSERT(r == NE || r == PE, "Cannot get partial electrode reaction current, only negative (NE) and positive electrodes (PE) are supported.");

   Vector amask({r == NE ? AN * LNE / NNE * NX : 0., 0., r == PE ? AP * LPE / NPE * NX : 0.});
   PWConstCoefficient a(amask);
   ProductCoefficient ajex(a, *_jex);
   TransformedCoefficient I(&ajex, _op, [=](real_t ajex, real_t op) { return ajex * exp( sign * 0.5 * op ); });
   QuadratureSpace x_qspace(_x_fespace->GetParMesh(), 2 * _x_fespace->FEColl()->GetOrder());
   return x_qspace.Integrate(I);
}

//
// Reaction Current for each particle
//

Array<real_t> EChemOperator::GetParticleReactionCurrent()
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

const real_t & EChemOperator::GetReactionCurrent(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get constant reaction current, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot get constant reaction current, only negative (NE) and positive electrodes (PE) are supported.");
   return _j->Eval()(r);
}

void EChemOperator::ConstructReactionCurrent()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         _j = new ReactionCurrentCoefficient();
         break;
      case P2D:
         _j = new ReactionCurrentCoefficient(T, *_jex, *_op);
         break;
   }
}

//
// Exchange Current
//

const real_t & EChemOperator::GetExchangeCurrent(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get constant exchange current, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot get constant exchange current, only negative (NE) and positive electrodes (PE) are supported.");
   return _jex->Eval()(r);
}

void EChemOperator::ConstructExchangeCurrent()
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

const real_t & EChemOperator::GetOpenCircuitPotential(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get constant open circuit potential, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot get constant open circuit potential, only negative (NE) and positive electrodes (PE) are supported.");
   return _ocp->Eval()(r);
}

void EChemOperator::ConstructOpenCircuitPotential()
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

const real_t & EChemOperator::GetOverPotential(const Region &r)
{
   MFEM_ASSERT(M == SPM || M == SPMe, "Cannot get constant overpotential, only SPM and SPMe are supported.");
   MFEM_ASSERT(r == NE || r == PE, "Cannot get constant overpotential, only negative (NE) and positive electrodes (PE) are supported.");
   return _op->Eval()(r);
}

void EChemOperator::ConstructOverPotential()
{
   switch (M)
   {
      case SPM:
      case SPMe:
         _op = new OverPotentialCoefficient(T, *_jex);
         break;
      case P2D:
         _op = new OverPotentialCoefficient(GetReferencePotential(E), GetReferencePotential(PE), *_sp_gf, *_ep_gf, *_ocp);
         break;
   }
}

//
// Voltage
//

real_t EChemOperator::GetVoltage()
{
   real_t Ve = (M == SPMe) ? GetVoltageMarquisCorrection() : 0;

   // Definition from JuBat: https://doi.org/10.1016/j.est.2023.107512
   real_t V = (GetOpenCircuitPotential(PE) - GetOpenCircuitPotential(NE) +
               GetOverPotential(PE) - GetOverPotential(NE) - Ve) * phi_scale;

   if (Mpi::Root())
   {
      std::cout << "[Rank " << Mpi::WorldRank() << "]"
                << " Voltage = " << V << std::endl;

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
            << V
            << std::endl;
   }

   _sc[1]->DebuggingCheck(_x);

   return V;
}

real_t EChemOperator::GetVoltageMarquisCorrection()
{
   GridFunctionCoefficient ec_gfc(_ec_gf);
   PWCoefficient ce_pwc;

   QuadratureSpace x_qspace(_x_fespace->GetParMesh(), 2 * _x_fespace->FEColl()->GetOrder());

   ce_pwc.UpdateCoefficient(NE, ec_gfc);
   real_t ce_ne_int = x_qspace.Integrate(ce_pwc) * (NX / NNE);
   ce_pwc.ZeroCoefficient(NE);

   ce_pwc.UpdateCoefficient(PE, ec_gfc);
   real_t ce_pe_int = x_qspace.Integrate(ce_pwc) * (NX / NPE);
   ce_pwc.ZeroCoefficient(PE);

   real_t eta_c = (2.0 * T / CE0) * (1 - TPLUS) * (ce_ne_int - ce_pe_int);
   real_t dphie = (I / KS) * (LNE / BNE / 3.0  + LSEP / BSEP + LPE / BPE / 3.0);
   real_t dphis = I / 3 * (LNE / SIGN + LPE / SIGP);

   return eta_c + dphie + dphis;
}

void EChemOperator::GetParticleDofs(Array<int> & particle_dofs, Array<Region> & particle_regions, Array<int> & particle_offsets)
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
   if (!IPAR)
      std::copy(sep_gdofs_set.begin(), sep_gdofs_set.end(), sep_gdofs.begin());

   Array<int> all_sep_gdofs(max_sep_gdofs * Mpi::WorldSize());
   MPI_Allgather(sep_gdofs.GetData(), max_sep_gdofs, MPI_INT,
                 all_sep_gdofs.GetData(), max_sep_gdofs, MPI_INT, MPI_COMM_WORLD);

   Array<int> boundary_dofs;
   if (!IPAR)
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
