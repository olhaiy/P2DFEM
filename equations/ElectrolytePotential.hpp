#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class ElectrolytePotential : public Equation
{
public:
  ElectrolytePotential(ParFiniteElementSpace & f) : Equation(f)
  {
    f.GetEssentialTrueDofs(Array<int>({1, 0}), ess_tdof_list);
  }
  virtual void Update(const BlockVector &, const Coefficient &) {}
  virtual void
  Update(const BlockVector & x, const GridFunctionCoefficient & ec_gfc, const Coefficient & j);
};
