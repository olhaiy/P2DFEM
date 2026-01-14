

#include "mfem.hpp"
#include "equations/Equation.hpp"

using namespace mfem;


#pragma once

class ElectrolytePotential : public Equation
{
   public:
      ElectrolytePotential(ParFiniteElementSpace &f) : Equation(f)
      {
         f.GetEssentialTrueDofs(Array<int>({1, 0}), ess_tdof_list);
      }
      void Update(const BlockVector &x, const Coefficient &j, const real_t &dt);
};
