#include <string>

namespace settings
{
extern const bool SPM;  // Whether running SPM model
extern const bool SPMe; // Whether running SPMe model
extern const bool P2D;  // Whether running P2D model

extern const unsigned NNE;  // Number of elements in the Negative Electrode
extern const unsigned NSEP; // Number of elements in the Separator
extern const unsigned NPE;  // Number of elements in the Positive Electrode
extern const unsigned
    NX; // Number of elements in the X-dimension (i.e Electrolye) (Sum of the above three)
extern const unsigned NR; // Number of elements in the R-dimension (i.e Particle)

extern const unsigned NNEPAR; // Number of Negative Electrode PARticle
extern const unsigned NPEPAR; // Number of Positive Electrode PARticle
extern const unsigned NPAR;   // Total number of  PARticles (Sum of the above two)

extern const unsigned NMACROP; // Number of Macro Potential equations //later change to NMACROPEQS!
extern const unsigned NMACROC; // Number of Macro Concentration equations
extern const unsigned NMACRO;  // Total Number of Macro Equations (Sum of the above two)
extern const unsigned NEQS;    // Total Number of Equations in the System

void init_settings(std::string m, int order);
}
