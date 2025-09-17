#include "mfem.hpp"
using namespace mfem;

class RegionalCurrent
{
   private:
      ParFiniteElementSpace * &_x_fespace;
      OverPotentialCoefficient * _op;
      ExchangeCurrentCoefficient * _jex;
      const Region &_r;
      TransformedCoefficient _i_tc;
      PWCoefficient _i_pwc;

   public:
      /// P2D
      RegionalCurrent(
      ParFiniteElementSpace * &x_fespace,
      const real_t & a,
      OverPotentialCoefficient & op,
      ExchangeCurrentCoefficient & jex,
      const Region & r):
      _x_fespace(x_fespace),
      _op(&op),        
      _jex(&jex),
      _r(r),
      _i_tc(_op, _jex, [=](real_t op, real_t jex) { return a * jex * exp( -0.5 * op); })
      {}

      /// P2D
      real_t operator()() 
      {
        QuadratureSpace x_qspace(_x_fespace->GetParMesh(), 2 * _x_fespace->FEColl()->GetOrder());
        _i_pwc.UpdateCoefficient(_r, _i_tc);
        return x_qspace.Integrate(_i_pwc);
      }
};
