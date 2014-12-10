#include <map>
#include <list>
#include <cstddef>
#include <datatools/properties.h>

#ifndef __N3AnaBase__gamma_tracking_h
#define __N3AnaBase__gamma_tracking_h 1

namespace gt {

  //! Implementation of the gamma tracking and combinator
  /*!
   * 2011-2012 Hugon Christophe CENBG
   *
   * The gamma tracking object is a combinator of integers in fonction of
   * doubles.
   * The integers represents PMs or vertexes number, the doubles are probabilities calculations.\n
   * It can give the longest combinaison or the one which has the highest probability.
   * It need an external probability calculator (tof_tools for NAT++) which is
   * able to give the probability between two PMs, and then it calculate the
   * combinaisons with chi square probabilies with (NbPMs-1) degrees of freedom.
   *
   * The most standard usage is :\n

   gamma_tracking gt;
   gt.AddProb (PM1, PM2, prob1);
   gt.AddProb (PM2, PM3, prob2);
   gt.AddProb (PM1, PM3, prob3);
   [...]
   gt.process ();
   refcoll_t gamma_tracked_coll=gt.get_reflects (lowest_prob);

   *
   */
  class gamma_tracking
  {
  public:

    /// collection of integer, ref for "reference of PM"
    typedef std::list <int> ref_t;

    /// collection of integer iterator, ref for "reference of PM"
    typedef std::list <int>::iterator ref_it;

    /// const collection of integer iterator, ref for "reference of PM"
    typedef std::list <int>::const_iterator ref_const_it;

    /// 2D collection of integer, ref for "reference of PM"
    typedef std::list<std::list <int> > refcoll_t;

    /// 2D collection of integer iterator, ref for "reference of PM"
    typedef std::list<std::list <int> >::iterator refcoll_it;

    /// 2D const collection of integer iterator, ref for "reference of PM"
    typedef std::list<std::list <int> >::const_iterator refcoll_const_it;

  public:

    /// Default constructor
    gamma_tracking();

    /// Copy constructor
    gamma_tracking(const gamma_tracking&);

    /// Default destructor
    ~gamma_tracking();

    /// Set the logging priority threshold
    void set_logging_priority(datatools::logger::priority priority_);

    /// Return the logging priority threshold
    datatools::logger::priority get_logging_priority() const;

    /// Check the initialization flag
    bool is_initialized() const;

    /// Set the initialization flag
    void set_initialized(bool);

    /// Initialization from parameters
    void initialize(const datatools::properties & config_);

    /// Check if contains calculated tracks
    bool has_tracks();

    /// Just add standalone gamma
    /*!< For complete gamma tracking, gamma_tracking::get_reflect gives also alone gammas*/
    void add(int number1_);

    /// Add 2 ref number combinaison with a proba in [0,1]
    void add_prob(int number1_, int number2_, double proba_);

    /// Add 2 ref number combinaison with a proba in [0,inf]
    void add_chi2(int number1_, int number2_, double chi2_);

    /// Add prestart.
    /*!< It impose starts during combinaison. Faster calculations,
      but less reliable than postarts. \sa gamma_tracking::get_reflects*/
    void add_start(int number_);

    /// cout information about current gamma tracking (old)
    void print();

    /// Counter of the current gamma tracking (old)
    void count();

    /// check if an element of gamma_tracking::ref_t values_ is in gamma_tracking::ref_t check_
    bool is_inside(const std::list <int>& check_, const std::list <int>& values_);

    /// check if value_ is in gamma_tracking::ref_t check_
    bool is_inside(const ref_t& check_, int value_);

    /// erase elements of gamma_tracking::ref_t values_ is in gamma_tracking::ref_t check_
    void extract(std::list <int>& source_, const std::list <int>& values_);

    /// to_ become a unique elements ref_t of from_+to_
    void put_inside(const std::list<int> &from_, std::list<int> &to_);

    /// if true, only by prob, else by size and then by prob.
    void set_absolute(bool a_);

    /// check gamma_tracking::set_absolute
    bool is_absolute();

    /// if true, forbid the starts to be elsewhere than at start
    void set_extern(bool e_);

    /// check gamma_tracking::set_extern
    bool is_extern();

    /// Set the minimal probability to continue next combinaisons
    void set_prob_min(double min_prob_);

    /*!<
      \param starts_ is the post start ref_t. After calculation, the function
      will return only gamma tracked starting with starts_.
      \param exclude_ is the ref_t to exclude from the calculations, and so, from the final result.
      \param deathless_starts_ means that the starts can be used in all gt.
      Usefull if starts represents vertexes instead of PM numbers.

      \attribute gamma_tracking::_extern_ allow only gamma tracked with starts, no elswhere\n
      gamma_tracking::_starts_ will be merged with starts_, and will be taken in care with extern and other

      \return The combinaison of gamma tracked. Each ref will be taken once
      (except for deathless starts), and the longest or the most probable gammas
      tracked will be given in gamma_tracking::ref_col_t format

      \sa gamma_tracking::set_absolute \sa gamma_tracking::set_extern
      \sa gamma_tracking::AddStart \sa gamma_tracking::process

    */
    /// Return the results
    refcoll_t get_reflects(double prob_, const ref_t* starts_ = NULL,
                           const ref_t *exclude_ = NULL, bool deathless_starts_ = false);

    /// Return all calculated combination.\sa gamma_tracking::process
    const refcoll_t get_all(); //removed the &

    /// Get the proba between two ref
    double get_prob(int scin_id1_, int scin_id2_) const;

    /// Get the proba for a gamma tracked
    double get_prob(const ref_t &scin_ids_) const;

    /// Get the chi square between two ref
    double get_chi2(int scin_id1_, int scin_id2_);

    /// Get the chi square for a gamma tracked
    double get_chi2(ref_t &scin_ids_);

    /// Classify the two ref_t in order of size and proba. Depend on _absolute_
    static bool sort_reflect(ref_t &ref1_, ref_t &ref2_);

    /// Classify a 2D ref_t in order of size and proba. Depend on _absolute_
    void sort_prob(std::list<std::list <int> > &list_);

    /// Get the chi square limit of _min_prob_ depend on degree of freedom
    double get_chi_limit(unsigned int);

    /// Main calculation before the gamma_tracking::get_reflects
    void process() ;

    /*!< Calculate all of the possible combinaisons of gamma tracked in the
      limit of gamma_tracking::_min_prob_. If there is prestart
      gamma_tracking::_starts_, it does the calculation only for combinaisons which starts with _starts_.
      \sa gamma_tracking::get_reflects \sa gamma_tracking::AddStart*/

    /// Reset the gamma tracking
    void reset() ;

    /// Tool to calculate each factorial only once
    static long factorial(int x);

    /* /\* interface i_serializable *\/ */
    /* virtual const std::string & get_serial_tag () const{}     */
    /* /\* interface i_clear *\/ */
    /* virtual void clear (){} */
    /* /\* interface i_tree_dumpable *\/ */
    /* virtual void tree_dump (std::ostream & out_         = std::clog,  */
    /* 			const std::string & title_  = "", */
    /* 			const std::string & indent_ = "", */
    /* 			bool inherit_               = false) const{} */

  private:

    /// Logging priority threshold
    datatools::logger::priority _logging_priority_;

    /// Initialization flag
    bool _initialized_;

    /// Prefer probability rather than size of gamma tracked
    bool _absolute_;

    /// Impose starts in the gamma tracked, not elsewhere
    bool _extern_;

    /// Maximum size of a gamma tracked
    int _max_;

    /// minimal probability to continue the combinating
    double _min_prob_;

    /// start collections (pre and post)
    ref_t _starts_;

    /// The full gamma tracked combinaisons
    refcoll_t _serie_;

    /// dictionnary of chi squares : deg of freedom: size-1 VS the chi2
    std::map<int,double> _min_chi2_;

    /// Dictionnary of chi square based on gamma tracked pointer
    std::map<const std::list<int>* ,double> _chi2_;

    /// Dictionnary of probabilities based on gamma tracked pointer
    std::map<const std::list<int>* ,double> _proba_;

    /// Factorial kept
    static std::map<long, double> _fact_;



  private:
    /*! What is saved in the file :
     * \param gamma_tracking::_serie_

     * \param gamma_tracking::_proba_

     * \param gamma_tracking::_max_

     * \param gamma_tracking::_starts_

     * \param gamma_tracking::_chi2_

     * \param gamma_tracking::_absolute_

     * \param gamma_tracking::_extern_

     * \param gamma_tracking::_min_prob_

     * \param gamma_tracking::_min_chi2_

     */
  };
}


#endif
