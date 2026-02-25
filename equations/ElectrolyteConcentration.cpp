
#include "equations/ElectrolyteConcentration.hpp"

void ElectrolyteConcentration::Update(const BlockVector &x, const GridFunctionCoefficient &ec_gfc, const Coefficient &j)
{
   // Mass coefficient.
   Vector mass_vec({
      /* NE */  EPS_N /* length scaling */  * (LNE  / NNE  * NX),
      /* SEP */ EPS_S /* length scaling */  * (LSEP / NSEP * NX),
      /* PE */  EPS_P /* length scaling */  * (LPE  / NPE  * NX)});
   PWConstCoefficient mass_part(mass_vec);
   ConstantCoefficient t_scale(te_scale);
   ProductCoefficient mass(mass_part, t_scale);

   // Source term.
   Vector source_vec({
      /* NE */ (1 - TPLUS) * AN /* length scaling */ * (LNE / NNE * NX),
      /* SEP */ 0.,
      /* PE */ (1 - TPLUS) * AP /* length scaling */ * (LPE / NPE * NX)});

   PWConstCoefficient source_part(source_vec);
   ProductCoefficient source(source_part, const_cast<Coefficient&>(j));

   // Diffusion coefficient.
   Vector D_scale_vec({
      /* NE */  BNE  /* length scaling */ / (LNE  / NNE  * NX),
      /* SEP */ BSEP /* length scaling */ / (LSEP / NSEP * NX),
      /* PE */  BPE  /* length scaling */ / (LPE  / NPE  * NX)});

   TransformedCoefficient D_coeff(&const_cast<GridFunctionCoefficient&>(ec_gfc), DE);
   PWConstCoefficient D_scale_coeff(D_scale_vec);
   ProductCoefficient D(D_scale_coeff, D_coeff);

   delete M;
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(mass));
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(D));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));

   Kmat.Mult(x.GetBlock(EC), b);
   b.Neg();
   b += Qvec;
}
