#include "mfem.hpp"
#include "cells/LGM50.hpp"

using namespace mfem;
using namespace LGM50;

namespace constants
{
enum PotentialBlock : int
{
  EPP, // Electrolyte Potential
  SPP, // Solid Potential (Postive and Negative Electrode region)
};

enum ConcentrationBlock : int
{
  ECC, // Electrolyte Concentration
  SCC  // Solid Concentration
};

enum Block : int
{
  P = 0, // Potential
  C = 2  // Concentration
};

enum XBlock : int
{
  EP = P + EPP, // Electrolyte Potential
  SP = P + SPP, // Solid Potential
  EC = C + ECC, // Electrolyte Concentration
  SC = C + SCC  // Solid Concentration
};

enum Region : int
{
  E,   // Electrolyte
  NE,  // Negative Electrode
  SEP, // Separator
  PE,  // Positive Electrode
  UNKNOWN
};

// Dimensional constants
const real_t F = 96485.33289; // Faraday constant, C/mol
const real_t R = 8.314;       // Universal gas constant, J/(mol*K)
const real_t T_ref = 298.;    // Reference temperature, K

// Scalings
const real_t t0 = 1.0;         // Time scale.
const real_t r0 = 1e-6;        // Length scale (particle)
const real_t L = ln + ls + lp; // Length scale (cell)

const real_t LNE = ln / L;  // Length of Negative Electrode
const real_t LSEP = ls / L; // Length of Separator
const real_t LPE = lp / L;  // Length of Positive Electrode

const real_t a0 = 1.0 / r0;

const real_t tn = F * cnmax * cell_area * L / I_typ; // Negative particle time scale.
const real_t tp = F * cpmax * cell_area * L / I_typ; // Positive particle time scale.

const real_t te = F * ce0 * cell_area * L / I_typ; // Electrolyte "particle" time scale.

const real_t te_scale = te / t0;

const real_t Dn_scale =
    r0 * r0 / t0; // Negative particle diffusion coefficient scale.  Units of m^2/s.
const real_t Dp_scale =
    r0 * r0 / t0; // Positive particle diffusion coefficient scale.  Units of m^2/s.

const real_t De_scale = L * L / te;

// Transport efficiency (inverse MacMullin number). This is B(x) in Planella, and is absorbed into
// the definition of kappa_ne/kappa_pe/kappa_sp in JuBat.
const real_t BPE = pow(eps_p, brugg);
const real_t BNE = pow(eps_n, brugg);
const real_t BSEP = pow(eps_s, brugg);

const real_t j_scale = I_typ / a0 / L / cell_area;

const real_t kn_scale = j_scale / cnmax / sqrt(ce0);
const real_t kp_scale = j_scale / cpmax / sqrt(ce0);

const real_t phi_scale = T_ref * R / F; // potential scale

// For scaling between particle time scale and cell time scale.
// Required for scaling the flux j in the SolidConcentration equation.
const real_t tn_scale = tn / t0; // time scale of negative electrode
const real_t tp_scale = tp / t0; // time scale of positive electrode

const real_t ce_scale = ce0;

const real_t sig_scale = L * I_typ / (cell_area * phi_scale);   // Electrode conductivity scale.
const real_t kappa_scale = L * I_typ / (cell_area * phi_scale); // Electrolyte conductivity scale.

// Scaled parameters
const real_t DN = Dn / Dn_scale; // scaled diffusion coefficient of each Negative particle
const real_t DP = Dp / Dp_scale; // scaled diffusion coefficient of each Positive particle

const real_t AN = An / a0; // scaled surface Area of each Negative particle // later SAN!
const real_t AP = Ap / a0; // scaled surface Area of each Positive particle

const real_t KN = kn_dim / kn_scale; // scaled reaction rate of each Negative particle
const real_t KP = kp_dim / kp_scale; // scaled reaction rate of each Positive particle

const real_t CN0 = cn0 / cnmax; // scaled initial Concentration of Negative particle
const real_t CP0 = cp0 / cpmax; // scaled initial Concentration of Positive particle

const real_t RN = rn / r0; // scaled Radius of Negative particle
const real_t RP = rp / r0; // scaled Radius of Positive particle

const real_t SIGP = sig_p / sig_scale; // Scaled positive electrode conductivity.
const real_t SIGN = sig_n / sig_scale; // Scaled negative electrode conductivity.

// Extras to be properly defined later.
const real_t CE0 = ce0 / ce_scale; // scaled initial Concentration of Electrolyte
const real_t I = 1.;               // scaled external current
const real_t T = 1.0;              // scaled Temperature.

const real_t EPS_P = eps_p;
const real_t EPS_N = eps_n;
const real_t EPS_S = eps_s;

const real_t TPLUS = tplus;

inline const real_t
UN(real_t ce)
{
  return Un(ce) / phi_scale;
}
inline const real_t
UP(real_t ce)
{
  return Up(ce) / phi_scale;
}
inline const real_t
DE(real_t ce)
{
  return De(ce * ce_scale) / De_scale;
}
inline const real_t
Kappa(real_t x)
{
  return kappa(x * ce_scale) / kappa_scale;
}

const real_t KS = Kappa(CE0); // / kappa_scale; // Scaled electrolyte conductivity.
}
