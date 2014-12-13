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
      gt::gamma_tracking::solution_type gamma_tracks;
      gt.get_reflects(1e-5, gamma_tracks);
      DT_LOG_DEBUG(get_logging_priority(), "Number of gammas = " << gamma_tracks.size());
      gt.print();

      for (gt::gamma_tracking::solution_type::const_iterator
             it = gamma_tracks.begin(); it != gamma_tracks.end(); ++it) {
        const gt::gamma_tracking::list_type & tmp_list = *it;
        //to change -> it_2 = it.begin() ?
        snemo::datamodel::particle_track::handle_type hPT(new snemo::datamodel::particle_track);
        track_.add_particle(hPT);
        hPT.grab().set_charge(snemo::datamodel::particle_track::neutral);

        for (gt::gamma_tracking::list_type::const_iterator jt = tmp_list.begin();
             jt != tmp_list.end(); ++jt) {
          if (jt == tmp_list.begin()) {
            double particle_time = -1;
            double particle_sigma_time = -1;
            double track_length = -1;
            double particle_energy=-1;
            double particle_sigma_energy=-1;
            double particle_x_vertex=-1;
            double particle_y_vertex=-1;
            double particle_z_vertex=-1;

            const snemo::datamodel::particle_track_data::particle_collection_type &
              the_particles = track_.get_particles();

            for (snemo::datamodel::particle_track_data::particle_collection_type::const_iterator
                   iparticle = the_particles.begin();
                 iparticle != the_particles.end();
                 ++iparticle) {

              const snemo::datamodel::particle_track & a_particle = iparticle->get();

              if (! a_particle.has_associated_calorimeter_hits()) {
                DT_LOG_DEBUG(get_logging_priority(),
                             "Particle track is not associated to any calorimeter block !");
                continue;
              }

              const snemo::datamodel::calibrated_calorimeter_hit::collection_type &
                the_calorimeters = a_particle.get_associated_calorimeter_hits ();

              if (the_calorimeters.size() > 1) {
                DT_LOG_DEBUG(get_logging_priority(),
                             "The particle is associated to more than 1 calorimeter !");
                continue;
              }

              if (a_particle.get_charge() == snemo::datamodel::particle_track::negative ||
                  a_particle.get_charge() == snemo::datamodel::particle_track::positive ||
                  a_particle.get_charge() == snemo::datamodel::particle_track::undefined) {

                // Look first if trajectory pattern is an helix or not
                const snemo::datamodel::tracker_trajectory & a_trajectory = a_particle.get_trajectory();
                const snemo::datamodel::base_trajectory_pattern & a_track_pattern = a_trajectory.get_pattern();
                const std::string & a_pattern_id = a_track_pattern.get_pattern_id();

                particle_time = the_calorimeters.at(0).get().get_time();
                particle_sigma_time = the_calorimeters.at(0).get().get_sigma_time();

                DT_THROW_IF(particle_time < 0, std::logic_error,
                            "The associated calorimeter time is negative !");

                DT_THROW_IF(particle_sigma_time < 0, std::logic_error,
                            "The associated calorimeter sigma time is negative !");

                particle_energy = the_calorimeters.at(0).get().get_energy();
                particle_sigma_energy = the_calorimeters.at(0).get().get_sigma_energy();


                DT_THROW_IF(particle_energy < 0, std::logic_error,
                            "The associated calorimeter energy is negative !");
                DT_THROW_IF(particle_sigma_energy < 0, std::logic_error,
                            "The associated calorimeter sigma energy is negative !");

                if (a_pattern_id == snemo::datamodel::line_trajectory_pattern::pattern_id()) {
                  const snemo::datamodel::line_trajectory_pattern * ptr_line
                    = dynamic_cast<const snemo::datamodel::line_trajectory_pattern *>(&a_track_pattern);
                  const geomtools::line_3d & a_line = ptr_line->get_segment();
                  track_length = a_line.get_length();
                } else if (a_pattern_id == snemo::datamodel::helix_trajectory_pattern::pattern_id()) {
                  const snemo::datamodel::helix_trajectory_pattern * ptr_helix
                    = dynamic_cast<const snemo::datamodel::helix_trajectory_pattern *>(&a_track_pattern);
                  const geomtools::helix_3d & a_helix = ptr_helix->get_helix();
                  track_length = a_helix.get_length();
                }

                const snemo::datamodel::particle_track::vertex_collection_type & the_vertices = a_particle.get_vertices();
                geomtools::vector_3d particle_source_vertex_position;
                geomtools::invalidate(particle_source_vertex_position);

                for (snemo::datamodel::particle_track::vertex_collection_type::const_iterator
                       ivertex = the_vertices.begin();
                     ivertex != the_vertices.end(); ++ivertex) {
                  const geomtools::blur_spot & a_vertex = ivertex->get();
                  const datatools::properties & aux = a_vertex.get_auxiliaries();

                  if(!snemo::datamodel::particle_track::vertex_is_on_source_foil(a_vertex)) {
                    DT_LOG_DEBUG(get_logging_priority(),
                                 "Vertex " << a_vertex.get_position()
                                 << " is not on the source foil !");
                    continue;
                  } else {
                    DT_LOG_WARNING(get_logging_priority(),
                                   "Possible vertex on foil !");
                    particle_source_vertex_position = a_vertex.get_position();

                    DT_THROW_IF(!geomtools::is_valid (particle_source_vertex_position), std::logic_error,
                                "Vertex position on source not valid !");

                    particle_x_vertex = particle_source_vertex_position.x();
                    particle_y_vertex = particle_source_vertex_position.y();
                    particle_z_vertex = particle_source_vertex_position.z();
                    break;
                  }
                }
              }

              double t_g = the_gamma_calos.at(*jt).time;
              double sigma_t_g = the_gamma_calos.at(*jt).sigma_time;
              double E_g = the_gamma_calos.at(*jt).energy;
              // double sigma_E_g = the_gamma_calos.at(*jt).sigma_energy;

              double gamma_track_length = sqrt(pow(particle_x_vertex-the_gamma_calos.at(*jt).position.x(),2) +
                                               pow(particle_y_vertex-the_gamma_calos.at(*jt).position.y(),2) +
                                               pow(particle_z_vertex-the_gamma_calos.at(*jt).position.z(),2));

              double t_g_th=gt::tof_computing::get_t_th(E_g,0.0,gamma_track_length);
              double t_e_th=gt::tof_computing::get_t_th(particle_energy,0.511,track_length);

              double sigma_particle_t_th = t_e_th*0.511*0.511/(particle_energy*(particle_energy+0.511)*(particle_energy+2*0.511))*particle_sigma_energy;

              double dt_int = particle_time-t_g-(t_e_th-t_g_th);

              double chi2_int = dt_int*dt_int/(particle_sigma_time*particle_sigma_time + sigma_t_g*sigma_t_g + sigma_particle_t_th*sigma_particle_t_th);
              double P_int =gsl_cdf_chisq_Q(chi2_int, 1);


              if (P_int>0.04) {
                DT_LOG_WARNING(get_logging_priority(),
                               "******Proba int : "<<P_int<<" -> Adding foil vertex !");

                snemo::datamodel::particle_track::handle_spot hBSv(new geomtools::blur_spot);
                hPT.grab().grab_vertices().push_back(hBSv);
                geomtools::blur_spot & spot_v = hBSv.grab();
                spot_v.set_hit_id(0);
                spot_v.grab_auxiliaries().store(snemo::datamodel::particle_track::vertex_type_key(),
                                                snemo::datamodel::particle_track::vertex_on_source_foil_label());
                spot_v.set_blur_dimension(geomtools::blur_spot::dimension_three);

                geomtools::vector_3d vertex_position (particle_x_vertex,particle_y_vertex,particle_z_vertex);

                spot_v.set_position(vertex_position);

              }

            }
          } //end if first calo from gamma

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
