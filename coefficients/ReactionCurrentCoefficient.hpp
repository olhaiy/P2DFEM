#include "mfem.hpp"
using namespace mfem;

class ReactionCurrentCoefficient : public Coefficient
{
private:
  ExchangeCurrentCoefficient * _jex = nullptr;
  OverPotentialCoefficient * _op = nullptr;

  Vector _j_vec;
  PWConstCoefficient _j_pwcc;
  TransformedCoefficient _j_tc;

public:
  /// SPM(e)
  ReactionCurrentCoefficient()
    : _j_vec({+I / AN / LNE, 0., -I / AP / LPE}),
      _j_pwcc(_j_vec),
      _j_tc(&_j_pwcc, [](real_t j) { return j; })
  {
  }

  /// P2D
  ReactionCurrentCoefficient(const real_t & T,
                             ExchangeCurrentCoefficient & jex,
                             OverPotentialCoefficient & op)
    : _jex(&jex),
      _op(&op),
      _j_tc(_jex, _op, [=](real_t jex, real_t op) { return 2 * jex * sinh(.5 * op / T); })
  {
  }

  /// SPM(e)
  virtual PWConstCoefficient & Eval() { return _j_pwcc; }

  /// P2D (and any integrators)
  virtual real_t Eval(ElementTransformation & T, const IntegrationPoint & ip) override
  {
    return _j_tc.Eval(T, ip);
  }
};
