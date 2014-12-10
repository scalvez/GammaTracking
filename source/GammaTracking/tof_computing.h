#ifndef GT_TOF_COMPUTING_H
#define GT_TOF_COMPUTING_H 1

#include "event.h"

namespace gt {

  class tof_computing
  {

  public :
    //   static const double kC=30;    //!< speed of light in cm.ns-1
    //   static const double kMe=511;   //!< electron mass in keV

    static double beta(double energy_, double mass_);

    static double get_t_th(double energy_, double mass_, double track_length_);

    /// Return distance between two calorimeter hits
    static double get_track_length(const event::calorimeter_hit & hit1_, const event::calorimeter_hit & hit2_);

    static double get_dt(const event::calorimeter_hit & hit1_, const event::calorimeter_hit & hit2_);

    static double get_chi2(const event::calorimeter_hit & hit1_, const event::calorimeter_hit & hit2_);

    static double get_proba(double chi2_, size_t ndf_ = 1);
  };

}
#endif // GT_TOF_COMPUTING_H
