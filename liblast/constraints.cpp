// constraints.cpp
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

#include "constraints.h"
#include "stats.h"

void ChisqLastConstraint::generateIntSubsets(set<int>& myset, set<set<int> >&subsets) {
    set<int>::iterator vi;
    for(vi = myset.begin(); vi != myset.end(); ++vi) {
      set<int> subset = myset;
      subset.erase(*vi);
      subsets.insert(subset);
      if(subset.size() > 1) generateIntSubsets(subset, subsets);
    }
}

float ChisqLastConstraint::ChiSqTest(float x, float y, unsigned int n_active, unsigned int n_inactive) {

  float res=0.0;
  // backup
  /*
  unsigned int na_tmp = na;
  unsigned int ni_tmp = ni;
  unsigned int n_tmp = n;
  bool activating_tmp=activating;

  na=n_active;
  ni=n_inactive;
  n=na+ni;
  res = gsl_cdf_chisq_P(ChiSq(x,y,true),1); // Always decide activating in this mode
  if (!activating) res*=-1.0;

  // restore
  na=na_tmp;
  ni=ni_tmp;
  n=n_tmp;
  activating=activating_tmp;
  */

  return res;
}


float ChisqLastConstraint::ChiSq(int x_val, vector<int> y) {
        assert(y.size() == nr_acts.size()); // equal class amounts as integrity constraint.
        int integrity = 0;
//        cout << endl;
        each(y) { 
          integrity+=y[i]; 
          //cout << "'" << y[i] << "'" << endl; 
        }
//        cout << "'" << integrity << "' '" << x_val << "'" <<  endl;
        assert(integrity == x_val);         // equal occurrence amounts as integrity constraint.

        float impact = 0.0;
        map<float, unsigned int>::iterator it;
        int i=0;

        impact = x_val/(float)n;
        chisq=0.0;
        for (it=nr_acts.begin(); it!=nr_acts.end(); it++) {
          float ev = it->second * impact;
          if (ev > 0) chisq += (y[i]-ev-0.5)*(y[i]-ev-0.5)/ev;
          i++;
        }
        return(chisq);
}

float ChisqLastConstraint::ChiSq(float x, float y, bool decide_activating) {
  /*
  float ea = 0.0, ei = 0.0, impact = 0.0;
  impact = x/(float)n;
  ea = na * impact; 
  ei = ni * impact; 
  if (decide_activating) {
    (y>ea) ? activating=1 : activating=0;
  }
  if (ea>0 && ei>0) chisq = (y-ea-0.5)*(y-ea-0.5)/ea + (x-y-ei-0.5)*(x-y-ei-0.5)/ei;
  */
  return(chisq);
}



float KSLastConstraint::KSTest(vector<float> all, vector<float> feat) {
  bool activating_tmp=activating;
  float res = KS(all, feat, 1);
  if (!activating) res*=-1.0;
  activating=activating_tmp;
  return res;
}

float KSLastConstraint::KS(vector<float> all_activities, vector<float> feat_activities, bool decide_activating) {

  // Kolmogorov-Smirnov Test
  // numerical recipies in C pp 626, bbrc_extended version with better sensitivity at the ends

  sort(feat_activities.begin(),feat_activities.end());
  sort(all_activities.begin(), all_activities.end());

  if (decide_activating) {
    // compute and compare medians to determine activation property
    if (!feat_activities.size()) {
      cerr << "Missing feature activities!" << endl;
      exit(1);
    }
    float aasum, aamedian, aamean, aavar, aadev, aaskew, aakurt;
    float fasum, famedian, famean, favar, fadev, faskew, fakurt;
    computeStats(all_activities.begin(), all_activities.end(), aasum, aamedian, aamean, aavar, aadev, aaskew, aakurt);
    computeStats(feat_activities.begin(), feat_activities.end(), fasum, famedian, famean, favar, fadev, faskew, fakurt);
    if (famedian > aamedian) activating=1; else activating = 0;
  }

  unsigned int j1=0, j2=0;
  float d,d1,d2,d_1,d_2,dt1,dt2,en1,en2,en,fn1=0,fn2=0,alam;
  d1 = d2 = d_1 = d_2 = 0.0;

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
