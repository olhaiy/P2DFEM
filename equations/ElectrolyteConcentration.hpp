#pragma once

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;

class ElectrolyteConcentration : public Equation
{
public:
  using Equation::Equation;
  virtual void Update(const BlockVector &, const Coefficient &) {}
  virtual void
  Update(const BlockVector & x, const GridFunctionCoefficient & ec_gfc, const Coefficient & j);
};
