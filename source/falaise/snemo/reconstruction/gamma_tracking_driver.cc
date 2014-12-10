/// \file falaise/snemo/reconstruction/gamma_tracking_driver.cc

// Ourselves:
#include <snemo/reconstruction/gamma_tracking_driver.h>

// Third party:
// - Boost:
#include <boost/fusion/iterator/next.hpp>
// - Bayeux/geomtools:
#include <geomtools/manager.h>

// This project:
#include <falaise/snemo/datamodels/particle_track.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/geometry/locator_plugin.h>
#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/xcalo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/datamodels/base_trajectory_pattern.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>
#include <falaise/snemo/datamodels/helix_trajectory_pattern.h>


// Gamma_Tracking library
#include <GammaTracking/event.h>
#include <GammaTracking/tof_computing.h>
#include <GammaTracking/gamma_tracking.h>

#include <gsl/gsl_cdf.h>

namespace snemo {

  namespace reconstruction {

    const std::string & gamma_tracking_driver::gamma_tracking_id()
    {
      static const std::string _id("gamma_tracking");
      return _id;
    }

    void gamma_tracking_driver::set_initialized(const bool initialized_)
    {
      _initialized_ = initialized_;
      return;
    }

    bool gamma_tracking_driver::is_initialized() const
    {
      return _initialized_;
    }

    void gamma_tracking_driver::set_logging_priority(const datatools::logger::priority priority_)
    {
      _logging_priority_ = priority_;
      return;
    }

    datatools::logger::priority gamma_tracking_driver::get_logging_priority() const
    {
      return _logging_priority_;
    }


    void gamma_tracking_driver::set_geometry_manager(const geomtools::manager & gmgr_)
    {
      DT_THROW_IF(is_initialized(), std::logic_error,
                  "Driver is already initialized !");
      _geometry_manager_ = &gmgr_;
      return;
    }

    const geomtools::manager & gamma_tracking_driver::get_geometry_manager() const
    {
      DT_THROW_IF(! has_geometry_manager(), std::logic_error,
                  "No geometry manager is setup !");
      return *_geometry_manager_;
    }

    bool gamma_tracking_driver::has_geometry_manager() const
    {
      return _geometry_manager_ != 0;
    }

    void gamma_tracking_driver::_set_defaults_()
    {
      _initialized_ = false;
      return;
    }

    // Constructor
    gamma_tracking_driver::gamma_tracking_driver()
    {
      _set_defaults_();
      return;
    }

    // Destructor
    gamma_tracking_driver::~gamma_tracking_driver()
    {
      if (is_initialized()) {
        reset();
      }
      return;
    }

    // Initialize the gamma tracker through configuration properties
    void gamma_tracking_driver::initialize(const datatools::properties & setup_)
    {
      // DT_THROW_IF (is_initialized(), std::logic_error, "Driver '" << get_id() << "' is already initialized !");
      DT_THROW_IF (is_initialized(), std::logic_error, "Driver 'GammaTracking' is already initialized !");

      DT_THROW_IF(! has_geometry_manager(), std::logic_error, "Missing geometry manager !");
      DT_THROW_IF(! get_geometry_manager().is_initialized(), std::logic_error,
                  "Geometry manager is not initialized !");

      // Logging priority
      datatools::logger::priority lp = datatools::logger::extract_logging_configuration (setup_);
      DT_THROW_IF(lp == datatools::logger::PRIO_UNDEFINED, std::logic_error,
                  "Invalid logging priority level for geometry manager !");
      set_logging_priority(lp);

      // Get geometry locator plugin
      const geomtools::manager & geo_mgr = get_geometry_manager();
      std::string locator_plugin_name;
      if (setup_.has_key("locator_plugin_name")) {
        locator_plugin_name = setup_.fetch_string("locator_plugin_name");
      } else {
        // If no locator plugin name is set, then search for the first one
        const geomtools::manager::plugins_dict_type & plugins = geo_mgr.get_plugins();
        for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin();
             ip != plugins.end();
             ip++) {
          const std::string & plugin_name = ip->first;
          if (geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(plugin_name)) {
            DT_LOG_DEBUG(get_logging_priority(), "Find locator plugin with name = " << plugin_name);
            locator_plugin_name = plugin_name;
            break;
          }
        }
      }
      // Access to a given plugin by name and type :
      DT_THROW_IF(! geo_mgr.has_plugin(locator_plugin_name) ||
                  ! geo_mgr.is_plugin_a<snemo::geometry::locator_plugin>(locator_plugin_name),
                  std::logic_error,
                  "Found no locator plugin named '" << locator_plugin_name << "'");
      _locator_plugin_ = &geo_mgr.get_plugin<snemo::geometry::locator_plugin>(locator_plugin_name);

      set_initialized(true);
      return;
    }

    // Reset the gamma tracker
    void gamma_tracking_driver::reset()
    {
      _set_defaults_();
      return;
    }

    // Main tracking method
    int gamma_tracking_driver::_process_algo(const snemo::datamodel::calibrated_calorimeter_hit::collection_type & hits_,
                                             snemo::datamodel::particle_track_data & track_)
    {
      DT_LOG_TRACE(get_logging_priority(), "Entering...");

      gt::event an_event;
      gt::event::calorimeter_collection_type & the_gamma_calos = an_event.grab_calorimeters();

      size_t count = 0;
      for (snemo::datamodel::calibrated_calorimeter_hit::collection_type::const_iterator
             icalo = hits_.begin(); icalo != hits_.end(); ++icalo, ++count) {
        const snemo::datamodel::calibrated_calorimeter_hit & a_calo_hit = icalo->get();
        {
          gt::event::calorimeter_hit dummy_hit;
          the_gamma_calos.insert(std::make_pair(count,dummy_hit));
        }

        gt::event::calorimeter_hit & new_calo_hit = the_gamma_calos[count];

        const geomtools::geom_id & a_gid = a_calo_hit.get_geom_id();
        const snemo::geometry::calo_locator & calo_locator   = _locator_plugin_->get_calo_locator();
        const snemo::geometry::xcalo_locator & xcalo_locator = _locator_plugin_->get_xcalo_locator();
        const snemo::geometry::gveto_locator & gveto_locator = _locator_plugin_->get_gveto_locator();
        if (calo_locator.is_calo_block_in_current_module(a_gid)) {
          calo_locator.get_block_position(a_gid, new_calo_hit.position);
          new_calo_hit.label = snemo::datamodel::particle_track::vertex_on_main_calorimeter_label();
        } else if (xcalo_locator.is_calo_block_in_current_module(a_gid)) {
          xcalo_locator.get_block_position(a_gid, new_calo_hit.position);
          new_calo_hit.label = snemo::datamodel::particle_track::vertex_on_x_calorimeter_label();
        } else if (gveto_locator.is_calo_block_in_current_module(a_gid)) {
          gveto_locator.get_block_position(a_gid, new_calo_hit.position);
          new_calo_hit.label = snemo::datamodel::particle_track::vertex_on_gamma_veto_label();
        } else {
          DT_THROW_IF(true, std::logic_error,
                      "Current geom id '" << a_gid << "' does not match any scintillator block !");
        }

        new_calo_hit.time         = a_calo_hit.get_time();
        new_calo_hit.sigma_time   = a_calo_hit.get_sigma_time();
        new_calo_hit.energy       = a_calo_hit.get_energy();
        new_calo_hit.sigma_energy = a_calo_hit.get_sigma_energy();
      }

      if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) {
        DT_LOG_DEBUG(get_logging_priority(), "Event dump: " << an_event);
      }

      gt::gamma_tracking gt;
      for (gt::event::calorimeter_collection_type::const_iterator
             icalo = the_gamma_calos.begin(); icalo != the_gamma_calos.end(); ++icalo) {
        for (gt::event::calorimeter_collection_type::const_iterator
               jcalo = boost::next(icalo); jcalo != the_gamma_calos.end(); ++jcalo) {
          gt::event::calorimeter_collection_type::const_iterator it1
            = icalo->second < jcalo->second ? icalo : jcalo;
          gt::event::calorimeter_collection_type::const_iterator it2
            = icalo->second < jcalo->second ? jcalo : icalo;
          const double tof_chi2 = gt::tof_computing::get_chi2(it1->second, it2->second);
          const double tof_prob = gt::tof_computing::get_internal_probability(tof_chi2);
          DT_LOG_DEBUG(get_logging_priority(), "X²(" << it1->first << "->"
                       << it2->first << ") = " << tof_chi2 << ", P = " << tof_prob);
          gt.add_prob(it1->first, it2->first, tof_prob);
        }
      }
      gt.process();

      // To be changed by returning list by reference
      const std::list<std::list<int> > gamma_tracked_coll = gt.get_reflects(1e-5);
      DT_LOG_DEBUG(get_logging_priority(), "Number of gammas = " << gamma_tracked_coll.size());

      // Maybe only take the best one otherwise you will get several gamma for
      // all the different lists
      for (std::list<std::list<int> >::const_iterator
             it = gamma_tracked_coll.begin(); it != gamma_tracked_coll.end(); ++it) {
        const std::list<int> & tmp_list = *it;
        //to change -> it_2 = it.begin() ?
        snemo::datamodel::particle_track::handle_type hPT(new snemo::datamodel::particle_track);
        track_.add_particle(hPT);
        hPT.grab().set_charge(snemo::datamodel::particle_track::neutral);

        for(std::list<int>::const_iterator jt = tmp_list.begin(); jt != tmp_list.end(); ++jt) {
          snemo::datamodel::particle_track::handle_spot hBS(new geomtools::blur_spot);
          hPT.grab().grab_vertices().push_back(hBS);
          geomtools::blur_spot & spot = hBS.grab();
          spot.set_hit_id(*jt);
          spot.grab_auxiliaries().store(snemo::datamodel::particle_track::vertex_type_key(),
                                        the_gamma_calos.at(*jt).label);
          spot.set_blur_dimension(geomtools::blur_spot::dimension_three);
          spot.set_position(the_gamma_calos.at(*jt).position);
        }
      }

      // build a vertex
      return 0;
    }  // end of namespace reconstruction

  }  // end of namespace snemo
}
/* OCD support */
#include <datatools/object_configuration_description.h>
DOCD_CLASS_IMPLEMENT_LOAD_BEGIN(snemo::reconstruction::gamma_tracking_driver, ocd_)
{
  ocd_.set_class_name("snemo::reconstruction::gamma_tracking_driver");
  ocd_.set_class_description("A driver class for the Gamma_Tracking algorithm");
  ocd_.set_class_library("Falaise_Gamma_Tracking");
  ocd_.set_class_documentation("The driver manager for the Gamma Tracking algorithms\n"
                               "/todo What does the manager do ?"
                               );


  ocd_.set_validation_support(true);
  ocd_.lock();
  return;
}
DOCD_CLASS_IMPLEMENT_LOAD_END() // Closing macro for implementation
DOCD_CLASS_SYSTEM_REGISTRATION(snemo::reconstruction::gamma_tracking_driver,
                               "snemo::reconstruction::gamma_tracking_driver")
