#include "mfem.hpp"
using namespace mfem;

class OpenCircuitPotentialCoefficient : public Coefficient
{
private:
  GridFunctionCoefficient * _surface_concentration_gfc = nullptr;

  const std::function<real_t(real_t)> _un;
  const std::function<real_t(real_t)> _up;

  const real_t * _scn = nullptr;
  const real_t * _scp = nullptr;

  TransformedCoefficient _ocp_ne_tc;
  TransformedCoefficient _ocp_pe_tc;

  PWConstCoefficient _ocp_pwcc;
  PWCoefficient _ocp_pwc;

public:
  /// SPM(e)
  OpenCircuitPotentialCoefficient(const std::function<real_t(real_t)> & un,
                                  const std::function<real_t(real_t)> & up,
                                  const real_t & scn,
                                  const real_t & scp)
    : _un(un),
      _up(up),
      _scn(&scn),
      _scp(&scp),
      _ocp_ne_tc(nullptr, _un),
      _ocp_pe_tc(nullptr, _up),
      _ocp_pwcc(3)
  {
  }

  /// P2D
  OpenCircuitPotentialCoefficient(const std::function<real_t(real_t)> & un,
                                  const std::function<real_t(real_t)> & up,
                                  GridFunctionCoefficient & sc)
    : _surface_concentration_gfc(&sc),
      _un(un),
      _up(up),
      _ocp_ne_tc(_surface_concentration_gfc, _un),
      _ocp_pe_tc(_surface_concentration_gfc, _up),
      _ocp_pwc(Array<int>({NE, PE}),
               Array<Coefficient *>({static_cast<Coefficient *>(&_ocp_ne_tc),
                                     static_cast<Coefficient *>(&_ocp_pe_tc)}))
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval()
  {
    _ocp_pwcc(NE) = _un(*_scn);
    _ocp_pwcc(PE) = _up(*_scp);
    return _ocp_pwcc;
  }

  /// P2D
  virtual real_t Eval(ElementTransformation & T, const IntegrationPoint & ip) override
  {
    return _ocp_pwc.Eval(T, ip);
  }
};
