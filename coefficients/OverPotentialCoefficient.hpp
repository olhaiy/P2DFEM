#include "mfem.hpp"
using namespace mfem;

class OverPotentialCoefficient : public Coefficient
{
private:
  GridFunctionCoefficient * _solid_potential_gfc = nullptr;
  GridFunctionCoefficient * _electrolyte_potential_gfc = nullptr;

  ExchangeCurrentCoefficient * _jex = nullptr;
  OpenCircuitPotentialCoefficient * _ocp = nullptr;

  FunctionCoefficient _rp_ne_fc;
  FunctionCoefficient _rp_pe_fc;
  PWCoefficient _rp_pwc;

  SumCoefficient _rel_se_sc;
  SumCoefficient _abs_se_sc;

  PWConstCoefficient _op_pwcc;
  SumCoefficient _op_sc;

public:
  /// SPM(e)
  OverPotentialCoefficient(const real_t & T, ExchangeCurrentCoefficient & jex)
    : _jex(&jex),
      _rp_ne_fc([](const Vector &) { return 0; }),
      _rp_pe_fc([](const Vector &) { return 0; }),
      _rel_se_sc(0, *_solid_potential_gfc),
      _abs_se_sc(0, *_solid_potential_gfc),
      _op_pwcc(3),
      _op_sc(0, *_solid_potential_gfc)
  {
  }

  /// P2D
  OverPotentialCoefficient(const real_t & rpe,
                           const real_t & rpp,
                           GridFunctionCoefficient & sp,
                           GridFunctionCoefficient & ep,
                           OpenCircuitPotentialCoefficient & ocp)
    : _solid_potential_gfc(&sp),
      _electrolyte_potential_gfc(&ep),
      _ocp(&ocp),
      _rp_ne_fc([&](const Vector &) { return 0 - rpe; }),
      _rp_pe_fc([&](const Vector &) { return rpp - rpe; }),
      _rp_pwc(Array<int>({NE, PE}),
              Array<Coefficient *>({static_cast<Coefficient *>(&_rp_ne_fc),
                                    static_cast<Coefficient *>(&_rp_pe_fc)})),
      _rel_se_sc(*_solid_potential_gfc, *_electrolyte_potential_gfc, 1, -1),
      _abs_se_sc(_rel_se_sc, _rp_pwc),
      _op_sc(_abs_se_sc, *_ocp, 1, -1)
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval()
  {
    _op_pwcc(NE) = 2 * T * asinh(+I / AN / LNE / 2.0 / _jex->Eval()(NE));
    _op_pwcc(PE) = 2 * T * asinh(-I / AP / LPE / 2.0 / _jex->Eval()(PE));
    return _op_pwcc;
  }

  /// P2D
  virtual real_t Eval(ElementTransformation & T, const IntegrationPoint & ip) override
  {
    return _op_sc.Eval(T, ip);
  }
};
