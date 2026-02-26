#include "mfem.hpp"

using namespace mfem;

namespace settings
{
bool SPM = false;
bool SPMe = false;
bool P2D = false;

unsigned NNE = 0;
unsigned NSEP = 0;
unsigned NPE = 0;
unsigned NX = 0;
unsigned NR = 10;

unsigned NNEPAR = 0;
unsigned NPEPAR = 0;
unsigned NPAR = 0;

unsigned NMACROP = 2;
unsigned NMACROC = 1;
unsigned NMACRO = NMACROP + NMACROC;
unsigned NEQS = 0;

void
init_settings(std::string m, int order)
{
  std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c) { return std::tolower(c); });

  if (m == "spm")
    SPM = true;
  else if (m == "spme")
    SPMe = true;
  else if (m == "p2d" || m == "dfn")
    P2D = true;
  else
    mfem_error("Unrecognised model.");

  if (SPM)
    NNE = NSEP = NPE = 0;
  else if (SPMe || P2D)
    NNE = NSEP = NPE = 10;
  NX = NNE + NSEP + NPE;

  if (SPM || SPMe)
    NNEPAR = NPEPAR = 1;
  else if (P2D)
  {
    NNEPAR = NNE * order + 1;
    NPEPAR = NPE * order + 1;
  }
  NPAR = NNEPAR + NPEPAR;
  NEQS = NMACRO + NPAR;
}
}
