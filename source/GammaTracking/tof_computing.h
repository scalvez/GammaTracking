#ifndef GT_TOF_COMPUTING_H
#define GT_TOF_COMPUTING_H 1

#include "event.h"

namespace gt {

  /// \brief Basic class to host Time-Of-Flight related functions
  class tof_computing
  {
  public :

    static double beta(double energy_, double mass_);

    static double get_t_th(double energy_, double mass_, double track_length_);

    /// Return distance between two calorimeter hits
    static double get_track_length(const event::calorimeter_hit & hit1_, const event::calorimeter_hit & hit2_);

    static double get_dt(const event::calorimeter_hit & hit1_, const event::calorimeter_hit & hit2_);

    /// Compute X² value between two calorimeter hits
    static double get_chi2(const event::calorimeter_hit & hit1_, const event::calorimeter_hit & hit2_);

    /// Return the internal probability given the X² value
    static double get_internal_probability(double chi2_, size_t ndf_ = 1);
  };

}
#endif // GT_TOF_COMPUTING_H
