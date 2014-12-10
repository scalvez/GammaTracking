#include "gamma_tracking.h"
#include <gsl/gsl_cdf.h>
#include <datatools/logger.h>
#include <algorithm>
#include <iostream>

using namespace std;

namespace gt {

  gamma_tracking::gamma_tracking()
  {
    _max_ = 0;
    _min_prob_ = 1e-5;
    _min_chi2_[1] = gsl_cdf_chisq_Qinv (_min_prob_,1);
    _absolute_ = false;
    _extern_ = false;
    // list<int> tamp;
    // _serie_.push_back (tamp);

  }

  gamma_tracking::gamma_tracking(const gamma_tracking& gt_)
  {
    _absolute_=gt_._absolute_;
    _extern_=gt_._extern_;
    _max_=gt_._max_;
    _min_prob_=gt_._min_prob_;
    _starts_=gt_._starts_;
    _serie_=gt_._serie_;
    _min_chi2_=gt_._min_chi2_;
    _fact_=gt_._fact_;

    for (map <const list<int>*,double>::const_iterator mit= gt_._chi2_.begin ();
	 mit!=gt_._chi2_.end();
	 mit++)
      {
	refcoll_it it = find (_serie_.begin(), _serie_.end(), *(mit->first));
	if (it!=_serie_.end())
	  _chi2_[&(*it)]=mit->second;
      }

  for (map<const list<int>* ,double>::const_iterator mit= gt_._proba_.begin ();
       mit!=gt_._proba_.end();
       mit++)
      {
	refcoll_it it = find (_serie_.begin(), _serie_.end(), *(mit->first));
	if (it!=_serie_.end())
	  _proba_[&(*it)]=mit->second;
      }

  }

  gamma_tracking::~gamma_tracking()
  {
    reset();
  }

  void gamma_tracking::set_logging_priority(datatools::logger::priority priority_)
  {
    _logging_priority_ = priority_;
    return;
  }

  datatools::logger::priority gamma_tracking::get_logging_priority() const
  {
    return _logging_priority_;
  }

  bool gamma_tracking::is_initialized() const
  {
    return _initialized_;
  }

  void gamma_tracking::set_initialized(bool initialized_)
  {
    _initialized_ = initialized_;
    return;
  }

  void gamma_tracking::initialize(const datatools::properties & config_)
  {
    DT_LOG_DEBUG(get_logging_priority(), "Entering...");
    DT_THROW_IF (is_initialized(), std::logic_error, "Already initialized !");


    if (config_.has_flag("absolute")) {
      _absolute_ = true; //or use method set_absolute ?
    }

    if (config_.has_flag("extern")) {
      _extern_ = true;
    }

    if (config_.has_key("max")) {
      _max_ = config_.fetch_integer("max");
    }

    if (config_.has_key("min_prob")) {
      _min_prob_ = config_.fetch_integer("min_prob");
    }

    set_initialized(true);
    DT_LOG_DEBUG(get_logging_priority(), "Exiting.");
    return;
  }

  bool gamma_tracking::has_tracks()
  {
    return _serie_.size();
  }

  void gamma_tracking::add(int number_)
  {
    list<int> tamp;
    tamp.push_back (number_);

    if (find (_serie_.begin (), _serie_.end (), tamp) == _serie_.end ())
      {
	_serie_.push_back (tamp);
	_proba_[&(_serie_.back ())] = 1;
	_chi2_[&(_serie_.back ())] = 0;
        ++_max_;
      }
  }

  void gamma_tracking::add_prob(int number1_, int number2_, double proba_)
  {
    double chi2;
    list<int> tamp;

    add (number1_);
    add (number2_);

    // std::cout<<"proba "<< number1_<<"-"<<number2_<<"  "<<proba_<<std::endl;
    if (number1_ == number2_ || proba_ < _min_prob_ )
      {
        // std::cout<<"proba too low  "<< std::endl;
      return;
      }

    tamp.push_back (number1_);
    tamp.push_back (number2_);

    if (find (_serie_.begin (), _serie_.end (), tamp) != _serie_.end ())
      return;

    _serie_.push_back (tamp);

    chi2 = gsl_cdf_chisq_Qinv (proba_,1);
    _chi2_[&(_serie_.back ())] = chi2;
    _proba_[&(_serie_.back ())] = proba_;

  }

  void gamma_tracking::add_chi2(int number1_, int number2_, double chi2_)
  {
    double proba = gsl_cdf_chisq_Q (chi2_,1);
    list<int> tamp;

    add (number1_);
    add (number2_);

    if (number1_ == number2_ || chi2_ > get_chi_limit (1) )
      return;

    tamp.push_back (number1_);
    tamp.push_back (number2_);

    if (find (_serie_.begin (), _serie_.end (), tamp) != _serie_.end ())
      return;

    _serie_.push_back (tamp);

    proba = gsl_cdf_chisq_Q (chi2_,1);
    _chi2_[&(_serie_.back ())] = chi2_;
    _proba_[&(_serie_.back ())] = proba;
  }

  void gamma_tracking::add_start(int number_)
  {
    bool its_inside = false;
    for (refcoll_it rit = _serie_.begin () ; rit != _serie_.end () && !its_inside ; rit++)
      its_inside = is_inside ((*rit), number_);

    if (!its_inside)
      {
	clog << "gamma_tracking:WARNING: PM " << number_ << " is not in the collection. You should add it before. Ignored\n";
	return;
      }

    _starts_.push_back (number_);

  }

  void gamma_tracking::print()
  {
    list<list<int> >::iterator l_it=_serie_.begin ();
    while ( _serie_.end () != l_it){
      cout << "for list ";
      for (list<int>::iterator it=(*l_it).begin (); it!= (*l_it).end (); it++)
	cout << (*it) << ' ';
      cout << "proba is " << _proba_[&(*l_it)] <<endl;
      l_it++;
    }
  }

  void gamma_tracking::count()
  {
    map<int,int> the_size;
    list<list<int> >::iterator l_it = _serie_.begin ();
    while ( _serie_.end () != l_it)
      {
        if (!the_size.count ((*l_it).size ()))
          the_size[(*l_it).size ()] = 0;

        the_size[(*l_it).size ()]++;
        l_it++;
      }

    cout << "the_size :\n"<<endl;

    for (map <int, int>::iterator it = the_size.begin (); it != the_size.end (); it++)
      {
        cout << it->first << ' ' << it->second << " instead of " << (factorial (_max_)/factorial (_max_-it->first))
             << " with a max of " << _max_ << endl;
      }
  }

  bool gamma_tracking::is_inside(const list <int> &check_, int value_)
  {
    if (find (check_.begin (), check_.end (), value_) != check_.end ())
      return true;
    return false;
  }

  bool gamma_tracking::is_inside(const list <int> &check_, const list <int> &values_)
  {
    for (list<int>::const_iterator it = values_.begin (); it != values_.end (); it++)
      if (find (check_.begin (), check_.end (), (*it)) != check_.end ())
	return true;
    return false;
  }

  void gamma_tracking::extract(std::list <int>& source_, const std::list <int>& values_)
  {
    for (ref_const_it it = values_.begin (); it != values_.end (); it++)
      {
	ref_it to_extract = find (source_.begin (), source_.end (), *it);
	if (to_extract != source_.end ())
	  source_.erase (to_extract);
      }
  }

  void gamma_tracking::put_inside(const list <int> &from_, list <int> &to_)
  {
    for (ref_const_it it = from_.begin (); it != from_.end (); it++)
      to_.push_back (*it);
    to_.sort ();
    to_.unique ();
  }

  void gamma_tracking::set_absolute(bool a_)
  {
    _absolute_=a_;
  }

  bool gamma_tracking::is_absolute()
  {
    return _absolute_;
  }

  void gamma_tracking::set_extern(bool e_)
  {
    _extern_=e_;
  }

  bool gamma_tracking::is_extern()
  {
    return _extern_;
  }

  void gamma_tracking::set_prob_min(double min_prob_)
  {
    _min_prob_ = min_prob_;
    for(map<int,double>::iterator it = _min_chi2_.begin (); it != _min_chi2_.end (); it++)
      {
	it->second = gsl_cdf_chisq_Qinv (_min_prob_, it->first);
      }
  }

  gamma_tracking::refcoll_t gamma_tracking::get_reflects(double prob_, const ref_t *starts_, const ref_t *exclude_, bool deathless_starts_)
  {

    // std::cout<<"Entering...";
    refcoll_t to_return;
    ref_t to_exclude;

    if (exclude_)
      to_exclude=*exclude_;

    ref_t starts;
    if (starts_)
      starts = *starts_;

    put_inside (_starts_, starts);

    sort_prob (_serie_);

    for (refcoll_it rit = _serie_.begin (); rit != _serie_.end (); rit++)
      {
	if ((*rit).size () > _max_ - to_exclude.size ()
	    || is_inside ((*rit), to_exclude)
	    || prob_ > _proba_[&(*rit)])
	  continue;

	if ((! starts.size () || is_inside (starts, *((*rit).begin ()))) )
	  {
	    if (_extern_
		&& (!is_inside (starts, (*rit).front ()) ||
		    ((*rit).size () == 2 && is_inside (starts, (*rit).back ()))))
	      continue;

	    to_return.push_back (*rit);

	    if ( starts.size () && deathless_starts_)
	      {
		ref_t serie = (*rit);
		serie.pop_front ();
		put_inside (serie, to_exclude);
	      }
	    else
	      {
		put_inside ((*rit), to_exclude);
	      }
	  }
	else
	  {
	    if (!_extern_)
	      {
		put_inside ((*rit), to_exclude);
		if (starts.size () && deathless_starts_)
		  extract (to_exclude, starts);
	      }
	  }
      }

    /*for (ref_it ex_it=to_exclude.begin(); ex_it!=to_exclude.end(); ex_it++)
      cout << *ex_it << ' ';
    cout << endl;
    */
    return to_return;

  }

  const gamma_tracking::refcoll_t gamma_tracking::get_all()
  {
    return _serie_;
  }

  double gamma_tracking::get_prob(int scin_id1_, int scin_id2_ ) const
  {
    ref_t l1;
    l1.push_back (scin_id1_);
    l1.push_back (scin_id2_);
    refcoll_const_it it = find (_serie_.begin (), _serie_.end (), l1);
    if (it != _serie_.end ())
      {
	double to_return= _proba_.at(&(*it));
	return to_return;
      }
    return -1.;
  }

  double gamma_tracking::get_prob(const ref_t &scin_ids_) const
  {
    refcoll_const_it it = find (_serie_.begin (), _serie_.end (), scin_ids_);
    if (it != _serie_.end ())
      {
	double to_return= _proba_.at(&(*it));
	return to_return;
      }
    return -1.;
  }

  double gamma_tracking::get_chi2(int scin_id1_, int scin_id2_ )
  {
    ref_t l1;
    l1.push_back (scin_id1_);
    l1.push_back (scin_id2_);
    refcoll_it it = find (_serie_.begin (), _serie_.end (), l1);
    if (it != _serie_.end ())
      return _chi2_[&(*it)];
    return -1.;
  }

  double gamma_tracking::get_chi2(ref_t &scin_ids_)
  {
    refcoll_it it = find (_serie_.begin (), _serie_.end (), scin_ids_);
    if (it != _serie_.end ())
      return _chi2_[&(*it)];
    return -1.;
  }

  bool gamma_tracking::sort_reflect(ref_t &ref1_, ref_t &ref2_)
  {
    if (ref1_.size () <= ref2_.size () )
      return false;
    return true;
  }

  void gamma_tracking::sort_prob(list<list <int> > &list_)
  {
    if (list_.size () <= 1)
      return;

    bool has_changed=true;
    while (has_changed)
      {
        has_changed=false;
        list<list <int> >::iterator it1=list_.begin ();
        list<list <int> >::iterator it2=list_.begin ();
        it2++;
        while (it2!=list_.end() && !has_changed)
          {
            if ((*it1).size () > 1 && (*it2).size ()>1 &&
                _proba_[&(*it1)] < _proba_[&(*it2)] &&
                (_absolute_ || (*it1).size () <= (*it2).size ()))
              {
                has_changed=true;
                list_.splice (it1, list_, it2);
              }
            else
              {
                it1++;
                it2++;
              }
          }
      }
    return;
  }

  double gamma_tracking::get_chi_limit(unsigned int freedom_)
  {
    if (!(_min_chi2_.count(freedom_)))
      {
	_min_chi2_[freedom_]=gsl_cdf_chisq_Qinv (_min_prob_,freedom_);
      }

    return _min_chi2_[freedom_];
  }

  void gamma_tracking::process()
  {
    list<int> tamp_list;
    bool has_next = false;
    unsigned int first_loop = 1;

    list<list<int> >::iterator l_it1 = _serie_.begin ();
    while (l_it1 != _serie_.end ())
      {
        for (list<list<int> >::iterator l_it2 = _serie_.begin(); l_it2 != _serie_.end(); l_it2++)
          {
            if ( (*l_it2).size () == 2
                 && (*l_it1).size () > first_loop
                 && (*l_it1).back () == (*l_it2).front ()
                 && !(is_inside ((*l_it1), (*l_it2).back ()))
                 && (!_starts_.size () || is_inside (_starts_, (*(*l_it1).begin ()))) )
              {
                bool starts_in = false;

                if (_extern_)
                  {
                    for (list<int>::iterator it = _starts_.begin ();
                         it != _starts_.end ();
                         it++)
                      {
                        if (find (++((*l_it1).begin ()), (*l_it1).end (), *it) != (*l_it1).end ())
                          {
                            starts_in = true;
                            break;
                          }

                        if ((*l_it2).back () == *it)
                          {
                            starts_in = true;
                            break;
                          }
                      }
                  }

                if (starts_in) continue;

                int freedom = (*l_it1).size () + (*l_it2).size () - 2;
                double chi2 = _chi2_[&(*l_it1)] + _chi2_[&(*l_it2)];
                if ( chi2 < get_chi_limit (freedom) )
                  {
                    double the_prob = gsl_cdf_chisq_Q (chi2, freedom);

                    tamp_list = (*l_it1);
                    tamp_list.insert (tamp_list.end (), ++((*l_it2).begin ()), (*l_it2).end () );
                    if (find (_serie_.begin (), _serie_.end (), tamp_list) == _serie_.end ())
                      {
                        has_next = true;
                        _serie_.push_front (tamp_list);
                      }

                    tamp_list.clear ();
                    _proba_[&(_serie_.front ())] = the_prob; //_proba_[&(*l_it1)]*_proba_[&(*l_it2)];
                    _chi2_[&(_serie_.front ())] = chi2;
                  }
              }
          }

        l_it1++;

        if (l_it1 == _serie_.end () && has_next)
          {
            l_it1=_serie_.begin ();
            has_next=false;
            first_loop=2;
          }
      }

    _serie_.sort(sort_reflect);
  }

  void gamma_tracking::reset()
  {
    _serie_.clear ();
    _proba_.clear ();
    _starts_.clear ();
  }

  long gamma_tracking::factorial(int x)
  {
    if (!_fact_.count (x))
      {
        long fac = 1;
        for (int i = 2; i <= x; i++)
          fac *= i;
        _fact_[x] = fac;
      }

    return (long) _fact_[x];
  }

  map<long, double> gamma_tracking::_fact_;
}
