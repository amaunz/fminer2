// constraints.h
// (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010

/*
   This file is part of LibLast (liblast).

   LibLast is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   LibLast is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LibLast.  If not, see <http://www.gnu.org/licenses/>.
   */


#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <set>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include "legoccurrence.h"
#include "database.h"

namespace fm {
  extern LastDatabase* last_database;
}

class LastConstraint {};

class ChisqLastConstraint : public LastConstraint {
  public:

    unsigned int na, ni, n;          // counters for a, i, n
    unsigned int fa, fi;             // counters for occurrences of spec feature                      (cand. for making private)
    float sig;                       // significance threshold, not accessed (r,rw) inside this class (only constructor)
    float chisq;                     // chisq is test results, written on every test                  (cand. for making private)
    float p, u;                      // p, u are test results, written on every test                  (cand. for making ro)
    bool active;                     // whether test is active, not accessed (r,rw) inside this class (only constructor)

    set<LastTid> fa_set, fi_set;     //                                                               (cand. for making private)
    bool activating;                 // defaults to deactivating (0)                                  (cand. for making ro)

    ChisqLastConstraint (float sig) : na(0), ni(0), n(0), fa(0), fi(0), sig(sig), chisq(0.0), p(0.0), u(0.0), active(0), activating(0) {}

    //!< Calculate chi^2 of current and upper bound for chi^2 of more specific features (see Morishita and Sese, 2000)
    template <typename OccurrenceType>
      void Calc(vector<OccurrenceType>& legocc) {

        chisq = 0.0; p = 0.0; u = 0.0; // init chisq, p, u

        LastLegActivityOccurrence(legocc);
        fa = fa_set.size(); // fa is y(I) in Morishita and Sese
        fi = fi_set.size(); // fi is x(I)-y(I)  in Morishita and Sese

        // chisq_p for current feature
        p = ChiSq(fa+fi,fa,1);         // set p to chisq test result => P IS NOT A P-VALUE here (legacy).

        // upper bound u for chisq_p of more specific features
        float u1 = 0.0, u2 = 0.0;      
        u1 = ChiSq(fa,fa,0);                                  // upper bound at
        u2 = ChiSq(fi,0,0);                                   // max{ chisq (y(I), y(I)) ,
        u = u1; if (u2>u1) u = u2;                            //      chisq (x(I)-y(I),0) }

      }
    float ChiSqTest(float x, float y, unsigned int n_active, unsigned int n_inactive); // on-the-fly test

  private:

    //!< Calculates chi^2 value
    float ChiSq(float x, float y, bool decide_activating);

    //!< Counts occurrences of legs in active and inactive compounds
    template <typename OccurrenceType>
      void LastLegActivityOccurrence(vector<OccurrenceType>& legocc) {

        fa_set.clear();
        fi_set.clear();

        each (legocc) { 

          if (fm::last_database->trees[legocc[i].tid]->activity == 1) {
            fa_set.insert(fm::last_database->trees[legocc[i].tid]->orig_tid); 
          }

          else if (fm::last_database->trees[legocc[i].tid]->activity == 0) {
            fi_set.insert(fm::last_database->trees[legocc[i].tid]->orig_tid); 
          }

        }

      }



};

class KSLastConstraint : public LastConstraint {
  public:
    vector<float> all;               // activity values (all)
    vector<float> feat;              // activity values (feature), written on every test              (cand. for making private)
    float sig;                       // significance threshold, not accessed (r,rw) inside this class (only constructor)
    float p;                         // p-value, written on every test                                (cand. for making ro)

    set<LastTid> fa_set, fi_set;     //                                                               (cand. for making private)
    bool activating;                 // defaults to deactivating (0)                                  (cand. for making ro)

    KSLastConstraint (float sig) : sig(sig), p(0.0), activating(0) {}

    template <typename OccurrenceType>
    void Calc(vector<OccurrenceType>& legocc) { LastLegActivityOccurrence(legocc); p = KS(all,feat,1); }
    float KSTest(vector<float> all, vector<float> feat);

  private:
    float KS(vector<float> all_activities, vector<float> feat_activities, bool decide_activating);
    //!< Stores activities of occurrences of legs
    template <typename OccurrenceType>
      void LastLegActivityOccurrence(vector<OccurrenceType>& legocc) {
        fa_set.clear();
        fi_set.clear();

        feat.clear();
        each (legocc) {
          feat.push_back(fm::last_database->trees[legocc[i].tid]->activity);
          fa_set.insert(fm::last_database->trees[legocc[i].tid]->orig_tid); 
        }
      }

};



#endif
