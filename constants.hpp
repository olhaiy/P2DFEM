#include "mfem.hpp"

using namespace mfem;

namespace constants {
    enum PotentialBlock : int {
        EPP,                       // Electrolyte Potential
        SPP,                       // Solid Potential (Postive and Negative Electrode region)
    };

    enum ConcentrationBlock : int {
        ECC,                       // Electrolyte Concentration
        SCC                        // Solid Concentration
    };

    enum Block : int {
        P = 0,                    // Potential
        C = 2                     // Concentration
    };

    enum XBlock : int {
        EP = P + EPP,             // Electrolyte Potential
        SP = P + SPP,             // Solid Potential
        EC = C + ECC,             // Electrolyte Concentration
        SC = C + SCC              // Solid Concentration
    };
                     
    enum Region : int {
        E,                        // Electrolyte
        NE,                       // Negative Electrode
        SEP,                      // Separator
        PE,                       // Positive Electrode
        UNKNOWN
    };

    enum Model : int {
        SPM,
        SPMe,
        P2D
    };

    extern const Model M;

    extern const unsigned NNE;     // Number of elements in the Negative Electrode
    extern const unsigned NSEP;    // Number of elements in the Separator
    extern const unsigned NPE;     // Number of elements in the Positive Electrode
    extern const unsigned NX;      // Number of elements in the X-dimension (i.e Electrolye) (Sum of the above three)
    extern const unsigned NR;      // Number of elements in the R-dimension (i.e Particle)

    extern const unsigned NNEPAR;  // Number of Negative Electrode PARticle
    extern const unsigned NPEPAR;  // Number of Positive Electrode PARticle
    extern const unsigned NPAR;    // Total number of  PARticles (Sum of the above two)

    extern const unsigned NMACROP; // Number of Macro Potential equations //later change to NMACROPEQS!
    extern const unsigned NMACROC; // Number of Macro Concentration equations
    extern const unsigned NMACRO;  // Total Number of Macro Equations (Sum of the above two)
    extern const unsigned NEQS;    // Total Number of Equations in the System

    extern const bool IPAR;        // Whether we place Interface Particles at the electrode-collector/separator interfaces

    extern const real_t LSEP;      // Length of Separator
    extern const real_t LPE;       // Length of Positive Electrode
    extern const real_t LNE;       // Length of Negative Electrode

    extern const real_t AN;        // scaled surface Area of each Negative particle // later SAN!
    extern const real_t AP;        // scaled surface Area of each Positive particle

    extern const real_t DN;        // scaled Diffusion coefficient of each Negative particle
    extern const real_t DP;        // scaled Diffusion coefficient of each Positive particle
    extern real_t DE(real_t ce);

    extern const real_t KN;        // scaled reaction rate of each Negative particle
    extern const real_t KP;        // scaled reaction rate of each Positive particle
    extern const real_t KS;

    extern const real_t RN;        // scaled Radius of Negative particle
    extern const real_t RP;        // scaled Radius of Positive particle

    extern const real_t CN0;       // scaled initial Concentration of Negative particle
    extern const real_t CP0;       // scaled initial Concentration of Positive particle

    extern const real_t CE0;       // scaled initial Concentration of Electrolyte

    extern const real_t I;         // scaled external current
    extern const real_t T;         // scaled Temperature.

    extern const real_t phi_scale; // potential scale
    extern const real_t tn_scale;  // time scale of negative electrode
    extern const real_t tp_scale;  // time scale of positive electrode
    extern const real_t te_scale;

    extern const real_t ce_scale;

    extern const real_t TPLUS;
    extern real_t BPE;
    extern real_t BNE;
    extern real_t BSEP;

    extern real_t EPS_P, EPS_N, EPS_S;
    extern real_t SIGN, SIGP;

    real_t UN(real_t);
    real_t UP(real_t);

    real_t Kappa(real_t x);

    void init_params(std::string m, int order);
}
