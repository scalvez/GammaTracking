// Standard libraries
#include <algorithm>
#include <random>

// Third party:
// - Bayeux/datatools:
#include <datatools/properties.h>

// This project
#include "event.h"
#include "tof_computing.h"
#include "gamma_tracking.h"

std::random_device rd;
std::default_random_engine generator(rd());

void generate_calorimeters(event & event_)
{
  const size_t total_nbr_calos = 12;
  std::uniform_int_distribution<int> distribution(0,total_nbr_calos-1);

  const size_t nbr_calos = 5;
  const size_t nbr_gammas = 1;

  // Sanity check
  if (nbr_calos * nbr_gammas > total_nbr_calos) {
    std::cerr << "Too much gammas for the given number of calorimeters ("
              << total_nbr_calos << ")" << std::endl;
    return;
  }

  for (size_t ig = 0; ig < nbr_gammas; ++ig) {
    // Initial start time
    double time = 0.0;
    for (size_t ic = 0; ic < nbr_calos; ++ic) {
      event::calorimeter_collection_type & the_calos = event_.grab_calorimeters();

      while(true) {
        const size_t gid = distribution(generator);
        if (std::find_if(the_calos.cbegin(), the_calos.cend(),
                         [gid] (const event::calorimeter_hit & h_) { return gid == h_.id; })
            != the_calos.cend()) {
          continue;
        }
        {
          event::calorimeter_hit dummy;
          the_calos.push_back(dummy);
        }
        auto icalo = std::prev(the_calos.end());
        icalo->id = gid;

        const double angle = 2*M_PI/double(total_nbr_calos);
        const double radius = 100;
        icalo->x = radius*cos(gid*angle);
        icalo->y = radius*sin(gid*angle);
        icalo->z = 0.0;

        if (ic == 0) {
          std::uniform_real_distribution<double> urdt(0.0, 100.0);
          time += urdt(generator);
          std::uniform_real_distribution<double> urde(0.0, 1.0);
          icalo->energy = urde(generator);
        } else {
          time += tof_computing::get_track_length(*icalo, *std::prev(icalo))/30.;
          std::uniform_real_distribution<double> urde(0.0, std::prev(icalo)->energy);
          icalo->energy = urde(generator);
        }

        icalo->energy = 1.0;
        // Smearing in energy
        const double fwhm2sig = 1.0/(2*sqrt(2*log(2.0)));
        icalo->sigma_energy = 0.08 * fwhm2sig * sqrt(icalo->energy) ;
        // std::normal_distribution<double> nde(icalo->energy, icalo->sigma_energy);
        // icalo->energy = nde(generator);
        // Smearing in time
        icalo->sigma_time = 0.250/sqrt(icalo->energy); // ns
        std::normal_distribution<double> ndt(time, icalo->sigma_time);
        icalo->time = ndt(generator);
        break;
      }
    }
  }

  return;
}

int main()
{
  const size_t nbr_events = 1;

  for (size_t i_evt = 0; i_evt < nbr_events; i_evt++) {
    event a_event;
    generate_calorimeters(a_event);

    std::cout << "Event #" << i_evt << std::endl;
    std::cout << a_event << std::endl;

    event::calorimeter_collection_type & the_calos = a_event.grab_calorimeters();
    // Sort calo hits in time
    the_calos.sort();

    // Gamma tracking
    datatools::properties setup;
    // setup.store_boolean("use_absolute_path", true);
    gt::gamma_tracking gt;
    // gt.initialize(setup);

    gt.set_absolute(true);

    for (auto icalo = the_calos.cbegin(); icalo != the_calos.cend(); ++icalo) {
      for (auto jcalo = std::next(icalo); jcalo != the_calos.cend(); ++jcalo) {
        const double tof_chi2 = tof_computing::get_chi2(*icalo, *jcalo);
        //std::clog << "X²(" << icalo->id << "->" << jcalo->id << ") = " << tof_chi2 << std::endl;
        const double tof_prob = tof_computing::get_proba(tof_chi2, the_calos.size() - 1);
        std::clog << "P(" << icalo->id << "->" << jcalo->id << ")  = " << tof_prob << std::endl;
        gt.add_prob(icalo->id, jcalo->id, tof_prob);
      }
    }

    // gt.set_absolute(true);
    // Replace Combine par process
    gt.process();
    const std::list<std::list<int>> gamma_tracked_coll = gt.get_reflects(1e-5);
    gt.print();
    // gt.count();
    std::cout << "Number of gammas found = " << gamma_tracked_coll.size() << std::endl;
    for (auto icol : gamma_tracked_coll) {
      std::cout << "Size of collection : " << icol.size() <<std::endl;
      for (auto ilis : icol) {
        std::cout<< ilis << "->";
      }
      std::cout << std::endl;
    }
  }
  return 0;
}
