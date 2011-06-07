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

    ChisqBbrcConstraint (float sig) : n(0), sig(sig), chisq(0.0), p(0.0), u(0.0) {}

    template <typename OccurrenceType>
    void Calc(vector<OccurrenceType>& legocc) {

        chisq = 0.0; p = 0.0; u = 0.0;

        BbrcLegActivityOccurrence(legocc);
        vector<int> f_sizes;
        int f_sum=0; each(f_sizes) f_sum+=f_sizes[i];

        // chisq_p for current feature
        p = ChiSq(f_sum, f_sizes);

        // upper bound u for chisq_p of more specific features
        float u=0.0;
        each(f_sizes) {
          int remember = f_sizes[i]; f_sizes[i]=0;
          float current = Chisq(f_sum-remember,f_sizes);
          if (current > u) u=current;
          f_sizes[i]=remember;
        }
    
    }


    private:

    //!< Calculates chi^2 and upper bound values
    float ChiSq(float x, float y);

    //!< Counts occurrences of legs in active and inactive compounds
    template <typename OccurrenceType>
    void BbrcLegActivityOccurrence(vector<OccurrenceType>& legocc) {

      fa_set.clear();
      fi_set.clear();
      fa_map.clear();
      fi_map.clear();

      std::pair< set<BbrcTid>::iterator, bool > insert_ret;

      each (legocc) { 

        BbrcTid orig_tid = fm::bbrc_database->trees[legocc[i].tid]->orig_tid;

        if (fm::bbrc_database->trees[legocc[i].tid]->activity == 1) {
            fa_map.insert(make_pair(orig_tid,1)); // each occurrence with 1, failure if present
            insert_ret = fa_set.insert(orig_tid); 
            if (!insert_ret.second) fa_map[orig_tid]++; // increase if present
        }

        else if (fm::bbrc_database->trees[legocc[i].tid]->activity == 0) {
            fi_map.insert(make_pair(orig_tid,1)); // each occurrence with 1, failure if present
            insert_ret = fi_set.insert(orig_tid); 
            if (!insert_ret.second) fi_map[orig_tid]++; // increase if present
        }

      }

    }
    


};

class KSBbrcConstraint : public BbrcConstraint {
    public:
    vector<float> all;
    vector<float> feat;
    float sig, p;
    set<BbrcTid> fa_set, fi_set;
    map<BbrcTid,int> fa_map, fi_map; // Store number of occurrences

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

      fa_set.clear();
      fi_set.clear();
      fa_map.clear();
      fi_map.clear();
      feat.clear();

      std::pair< set<BbrcTid>::iterator, bool > insert_ret;
      each (legocc) {
        feat.push_back(fm::bbrc_database->trees[legocc[i].tid]->activity);
        BbrcTid orig_tid = fm::bbrc_database->trees[legocc[i].tid]->orig_tid;
        fa_map.insert(make_pair(orig_tid,1)); // each occurrence with 1, failure if present
        insert_ret = fa_set.insert(orig_tid); 
        if (!insert_ret.second) fa_map[orig_tid]++; // increase if present
      }
    }

};

#endif
