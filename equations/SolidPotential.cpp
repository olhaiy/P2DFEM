
#include "equations/SolidPotential.hpp"

void SolidPotential::Update(const BlockVector &x, const Coefficient &j, const real_t &dt)
{
   // Source term.
   Vector source_vec({
      /* NE */ AN /* length scaling */ * (LNE / NNE * NX),
      /* SEP */ 0.,
      /* PE */ AP /* length scaling */ * (LPE / NPE * NX)});

   PWConstCoefficient source_part(source_vec);
   ProductCoefficient source(source_part, const_cast<Coefficient&>(j));

   // Effective conductivity (does not account for electrode filler).
   Vector sigma_vec({
      /* NE */ (1 - EPS_N) * SIGN /* length scaling */ / (LNE / NNE * NX),
      /* SEP */ 0.,
      /* PE */ (1 - EPS_P) * SIGP /* length scaling */ / (LPE / NPE * NX)});

   PWConstCoefficient sigma(sigma_vec);

   delete K;
   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(sigma));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);

   delete Q;
   Q = new ParLinearForm(&fespace);
   Q->AddDomainIntegrator(new DomainLFIntegrator(source));
   Q->Assemble();
   Qvec = std::move(*(Q->ParallelAssemble()));
   Qvec.SetSubVector(ess_tdof_list, 0.0); // do we need this?

   Kmat.Mult(x.GetBlock(SP), b);
   b += Qvec;
   b *= -1./dt;
}
