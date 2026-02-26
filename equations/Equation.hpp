#pragma once

#include "mfem.hpp"
#include "parameters/constants.hpp"
#include "parameters/settings.hpp"

using namespace mfem;
using namespace constants;
using namespace settings;

class Equation
{
protected:
  ParFiniteElementSpace & fespace;

  Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

  ParBilinearForm * M = nullptr;
  ParBilinearForm * K = nullptr;
  ParLinearForm * Q = nullptr;

  HypreParMatrix Mmat;
  HypreParMatrix Kmat;
  HypreParVector * Qvec = nullptr;

  HypreSmoother prec; // Preconditioner for the implicit solver

  mutable Vector b; // auxiliary vector

public:
  Equation(ParFiniteElementSpace & f) : fespace(f), b(f.GetTrueVSize()) {};

  const HypreParMatrix & GetM() const { return Mmat; };
  const HypreParMatrix & GetK() const { return Kmat; };
  const Vector & GetZ() const { return b; };

  /// Update the diffusion BilinearForm K using the given true-dof vector `x`.
  virtual void Update(const BlockVector & x, const Coefficient & j) = 0;
  virtual void
  Update(const BlockVector & x, const GridFunctionCoefficient & u, const Coefficient & j) = 0;

  virtual ~Equation()
  {
    delete M;
    delete K;
    delete Q;
    delete Qvec;
  }
};
