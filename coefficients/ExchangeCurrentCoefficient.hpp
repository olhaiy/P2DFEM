#include "mfem.hpp"
using namespace mfem;

class ExchangeCurrentCoefficient : public Coefficient
{
private:
  GridFunctionCoefficient * _surface_concentration_gfc = nullptr;
  GridFunctionCoefficient * _electrolyte_concentration_gfc = nullptr;

  const real_t * _scn = nullptr;
  const real_t * _scp = nullptr;

  TransformedCoefficient _jex_ne_tc;
  TransformedCoefficient _jex_pe_tc;

  Vector _jex_vec;
  PWConstCoefficient _jex_pwcc;
  PWCoefficient _jex_pwc;

public:
  /// SPM
  ExchangeCurrentCoefficient(const real_t & kn,
                             const real_t & kp,
                             const real_t & scn,
                             const real_t & scp,
                             const real_t & ec)
    : _scn(&scn),
      _scp(&scp),
      _jex_ne_tc(nullptr, [](real_t) { return 0; }),
      _jex_pe_tc(nullptr, [](real_t) { return 0; }),
      _jex_vec({kn * sqrt(ec), 0., kp * sqrt(ec)}),
      _jex_pwcc(3)
  {
  }

  /// SPMe
  ExchangeCurrentCoefficient(const real_t & kn,
                             const real_t & kp,
                             const real_t & scn,
                             const real_t & scp,
                             GridFunctionCoefficient & ec)
    : _electrolyte_concentration_gfc(&ec),
      _scn(&scn),
      _scp(&scp),
      _jex_ne_tc(_electrolyte_concentration_gfc, [=](real_t ec) { return kn * sqrt(ec); }),
      _jex_pe_tc(_electrolyte_concentration_gfc, [=](real_t ec) { return kp * sqrt(ec); })
  {
  }

  /// P2D
  ExchangeCurrentCoefficient(const real_t & kn,
                             const real_t & kp,
                             GridFunctionCoefficient & sc,
                             GridFunctionCoefficient & ec)
    : _surface_concentration_gfc(&sc),
      _electrolyte_concentration_gfc(&ec),
      _jex_ne_tc(_surface_concentration_gfc,
                 _electrolyte_concentration_gfc,
                 [=](real_t sc, real_t ec) { return kn * sqrt(sc * ec * (1 - sc)); }),
      _jex_pe_tc(_surface_concentration_gfc,
                 _electrolyte_concentration_gfc,
                 [=](real_t sc, real_t ec) { return kp * sqrt(sc * ec * (1 - sc)); }),
      _jex_pwc(Array<int>({NE, PE}),
               Array<Coefficient *>({static_cast<Coefficient *>(&_jex_ne_tc),
                                     static_cast<Coefficient *>(&_jex_pe_tc)}))
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval()
  {
    /// SPMe
    if (_electrolyte_concentration_gfc)
    {
      ParFiniteElementSpace * x_h1space =
          static_cast<const ParGridFunction *>(_electrolyte_concentration_gfc->GetGridFunction())
              ->ParFESpace();
      QuadratureSpace x_qspace(x_h1space->GetParMesh(), 2 * x_h1space->FEColl()->GetOrder());

      /// NE
      _jex_pwc.UpdateCoefficient(NE, _jex_ne_tc);
      real_t integral_ne = x_qspace.Integrate(_jex_pwc);
      _jex_pwc.ZeroCoefficient(NE);

      /// PE
      _jex_pwc.UpdateCoefficient(PE, _jex_pe_tc);
      real_t integral_pe = x_qspace.Integrate(_jex_pwc);
      _jex_pwc.ZeroCoefficient(PE);

      Vector c({integral_ne * sqrt(*_scn * (1 - *_scn)) / NNE * NX,
                0.,
                integral_pe * sqrt(*_scp * (1 - *_scp)) / NPE * NX});
      _jex_pwcc.UpdateConstants(c);
    }
    /// SPM
    else
    {
      _jex_pwcc(NE) = _jex_vec(NE - 1) * sqrt(*_scn * (1 - *_scn));
      _jex_pwcc(PE) = _jex_vec(PE - 1) * sqrt(*_scp * (1 - *_scp));
    }

    return _jex_pwcc;
  }

  /// P2D
  virtual real_t Eval(ElementTransformation & T, const IntegrationPoint & ip) override
  {
    return _jex_pwc.Eval(T, ip);
  }
};
