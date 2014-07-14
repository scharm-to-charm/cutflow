#ifndef CTAG_DEFS_HH
#define CTAG_DEFS_HH

// enums used by the calibration tool
namespace ctag {
  enum Flavor { B, C, U, T, DATA};
}

// log(pc / pu) selection minimum
const double JFC_MEDIUM_ANTI_U_CUT =  0.95;

// log(pc / pb) selection minimum
const double JFC_MEDIUM_ANTI_B_CUT = -0.90;


#endif
