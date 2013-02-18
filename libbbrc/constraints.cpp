// constraints.cpp
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

#include "constraints.h"


void ChisqBbrcConstraint::generateIntSubsets(set<int>& myset, set<set<int> >&subsets) {
    set<int>::iterator vi;
    for(vi = myset.begin(); vi != myset.end(); ++vi) {
      set<int> subset = myset;
      subset.erase(*vi);
      subsets.insert(subset);
      if(subset.size() > 1) generateIntSubsets(subset, subsets);
    }
}

float ChisqBbrcConstraint::ChiSq(int x_val, vector<int> y) {
        assert(y.size() == nr_acts.size()); // equal class amounts as integrity constraint.
        int integrity = 0;
        each(y) integrity+=y[i];
        assert(integrity == x_val);         // equal occurrence amounts as integrity constraint.

        int i=0;
        float impact = 0.0;
        map<float, unsigned int>::iterator it;

        impact = x_val/(float)n;
        chisq=0.0;
        for (it=nr_acts.begin(); it!=nr_acts.end(); it++) {
          float ev = it->second * impact;
          if (ev > 0) chisq += (y[i]-ev-0.5)*(y[i]-ev-0.5)/ev;
          i++;
        }
        chisq=chisq/1000; // AM: Quick hack dividing chisq/1000 - remove when weights have been converted to floats
        return(chisq);
}

float KSBbrcConstraint::KS(vector<float> all_activities, vector<float> feat_activities) {

    // Kolmogorov-Smirnov Test
    // numerical recipies in C pp 626, bbrc_extended version with better sensitivity at the ends

    sort(feat_activities.begin(),feat_activities.end());
    sort(all_activities.begin(), all_activities.end());

    unsigned int j1,j2;
    j1 = j2 = 0;

    float d,d1,d2,d_1,d_2,dt1,dt2,en1,en2,en,fn1,fn2,alam;
    d1 = d2 = d_1 = d_2 = dt1 = dt2 = en1 = en2 = en = fn1 = fn2 = alam = 0.0;

    en1 = all_activities.size();
    en2 = feat_activities.size();

    while (j1 < en1 && j2 < en2) {
        if ((!isnan(d1=all_activities[j1])) && (!isnan(d2=feat_activities[j2]))) {
            if (d1 <= d2) { // next step is in all_activities
                j1++;
                fn1=j1/en1;
            }

            if (d2 <= d1) { // next step is in feat_activities
                j2++;
                fn2=j2/en2;
            }

        }
        else {
            if (isnan(d1)) {
                j1++;
            }
            if (isnan(d2)){
                j2++;
            }
            continue;
        }
        dt1=fn2-fn1;
        dt2=fn1-fn2;
        if (dt1 > d_1) d_1=dt1;
        if (dt2 > d_2) d_2=dt2;
    }
    d = d_1 + d_2;
    en=sqrt(en1*en2/(en1+en2));
    alam=(en+0.155+0.24/en)*d;


    // KS probability function
    float p,a2,fac=2,sum=0,term,termbf=0,fac2;
    a2 = -2.0 *alam*alam;
    p = 0;
    fac2 = 4.0*alam*alam;
    for (int j=1;j<=100;j++) {
        term = fac*((fac2*j*j)-1.0)*exp(a2*j*j);
        sum += term;
        if (fabs(term) <= 0.001*termbf || fabs(term) <= 1.0e-8*sum) {
            p=1-sum;
            break;
        }
        termbf=fabs(term);
    }
    return p;
}
