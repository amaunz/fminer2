// constraints.h
// (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010

/*
    This file is part of LibBbrc (libbbrc).

    LibBbrc is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LibBbrc is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LibBbrc.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <set>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include "legoccurrence.h"
#include "database.h"
#include <assert.h>

namespace fm {
    extern BbrcDatabase* bbrc_database;
}

class BbrcConstraint {};

class ChisqBbrcConstraint : public BbrcConstraint {
    public:
    map<float, unsigned int> nr_acts;
    unsigned int n;
    unsigned int fa, fi;
    float sig, chisq, p, u;
    bool active;
    map<float, set<BbrcTid> > f_sets;
    map<float, map<BbrcTid,int> > f_maps; 
    map<int, float> df_thresholds;

    ChisqBbrcConstraint (float sig) : n(0), sig(sig), chisq(0.0), p(0.0), u(0.0) {
      df_thresholds[1]=3.84;
      df_thresholds[2]=5.99;
      df_thresholds[3]=7.82;
      df_thresholds[4]=9.49;
    }

    template <typename OccurrenceType>
    void Calc(vector<OccurrenceType>& legocc) {

        chisq = 0.0; p = 0.0; u = 0.0;

        BbrcLegActivityOccurrence(legocc);
        vector<int> f_sizes;
        map<float, set<BbrcTid> >::iterator  f_sets_it;
        for (f_sets_it=f_sets.begin(); f_sets_it!=f_sets.end(); f_sets_it++) f_sizes.push_back(f_sets_it->second.size());
        int f_sum=0; each(f_sizes) f_sum+=f_sizes[i];
        // chisq_p for current feature
        p = ChiSq(f_sum, f_sizes);

        // upper bound u for chisq_p of more specific features
        u=0.0;
        set<int> f_indices;
        for (int i=0; i<f_sizes.size(); i++) f_indices.insert(i);
        set<set<int> > f_subsets;
        set<set<int> >::iterator f_subsets_it;
        generateIntSubsets(f_indices, f_subsets);
        for (f_subsets_it = f_subsets.begin(); f_subsets_it!=f_subsets.end(); f_subsets_it++) {
          if (f_subsets_it->size() > 0) {
            f_sum=0;
            vector<int> f_selected_sizes;
            for (int j=0; j<f_sizes.size(); j++) { 
              if (f_subsets_it->find(j)!=f_subsets_it->end()) {
                f_selected_sizes.push_back(f_sizes[j]);
                f_sum+=f_sizes[j];
              }
              else f_selected_sizes.push_back(0);
            }
            float current = ChiSq(f_sum,f_selected_sizes);
            if (current > u) u=current;
          }
       }
//       cout << "U: " << u << " P: " << p << endl;
    }


    private:

    //!< Calculates chi^2 and upper bound values
    float ChiSq(int x_val, vector<int> y);
    void generateIntSubsets(set<int>& myset, set<set<int> >&subsets);
 

    //!< Counts occurrences of legs in active and inactive compounds
    template <typename OccurrenceType>
    void BbrcLegActivityOccurrence(vector<OccurrenceType>& legocc) {

      f_sets.clear();
      f_maps.clear();

      std::pair< set<BbrcTid>::iterator, bool > insert_ret;

      map<float, unsigned int>::iterator nr_acts_it;
      for (nr_acts_it = nr_acts.begin(); nr_acts_it != nr_acts.end(); nr_acts_it++) {
        set<BbrcTid> tmp;
        f_sets[nr_acts_it->first]=tmp;
      }

      each (legocc) { 
        float activity = fm::bbrc_database->trees[legocc[i].tid]->activity;
        BbrcTid orig_tid = fm::bbrc_database->trees[legocc[i].tid]->orig_tid;

        f_maps[activity].insert(make_pair(orig_tid,1)); // each occurrence with 1, failure if present
        insert_ret = f_sets[activity].insert(orig_tid); 
        if (!insert_ret.second) f_maps[activity][orig_tid]++; // increase if present

      }

    }
    


};

class KSBbrcConstraint : public BbrcConstraint {
    public:
    vector<float> all;
    vector<float> feat;
    float sig, p;
    map<float, set<BbrcTid> > f_sets;
    map<float, map<BbrcTid,int> > f_maps; 

    KSBbrcConstraint (float sig) : sig(sig), p(0.0) {}

    template <typename OccurrenceType>
    void Calc(vector<OccurrenceType>& legocc) {
        BbrcLegActivityOccurrence(legocc);
        p = KS(all,feat);
    }

    private:
    //!< Calculates KS values
    float KS(vector<float> all_activities, vector<float> feat_activities);

    //!< Stores activities of occurrences of legs
    template <typename OccurrenceType>
    void BbrcLegActivityOccurrence(vector<OccurrenceType>& legocc) {

     feat.clear();
     f_sets.clear();
     f_maps.clear();

      std::pair< set<BbrcTid>::iterator, bool > insert_ret;
      each (legocc) {
        float activity = fm::bbrc_database->trees[legocc[i].tid]->activity;
        BbrcTid orig_tid = fm::bbrc_database->trees[legocc[i].tid]->orig_tid;

        f_maps[0.0].insert(make_pair(orig_tid,1)); // each occurrence with 1, failure if present (use 0.0 as dummy key for regression)
        insert_ret = f_sets[0.0].insert(orig_tid); 
        if (!insert_ret.second) { 
          f_maps[0.0][orig_tid]++; // increase if present
          feat.push_back(activity);
        }
      }
    }

};

#endif
