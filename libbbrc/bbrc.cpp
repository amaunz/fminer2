// bbrc.cpp
// modified 2010
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

#include "bbrc.h"
#include "globals.h"
#include "ranker.h"


// 1. Constructors and Initializers

Bbrc::Bbrc() : init_mining_done(false) {
  if (!fm::bbrc_instance_present) {
      fm::bbrc_database = NULL; fm::bbrc_statistics = NULL; fm::bbrc_chisq = NULL; fm::bbrc_result = NULL;
      Reset();
      Defaults();
      fm::bbrc_instance_present=true;
      if (getenv("FMINER_LAZAR")) fm::bbrc_do_yaml = false;
      if (getenv("FMINER_SMARTS")) fm::bbrc_gsp_out = false; 
      if (getenv("FMINER_PVALUES")) fm::bbrc_pvalues = true;
      if (getenv("FMINER_NO_AROMATIC_WC")) fm::bbrc_aromatic_wc = false;
      if (getenv("FMINER_SILENT")) {
        FILE* fp = freopen ("fminer_debug.txt","w",stderr);
      }
      if (getenv("FMINER_NR_HITS")) fm::bbrc_nr_hits = true;
  }
  else {
    cerr << "Error! Cannot create more than 1 instance." << endl; 
    exit(1);
  }
}

Bbrc::Bbrc(int _type, unsigned int _minfreq) : init_mining_done(false) {
  if (!fm::bbrc_instance_present) {
      fm::bbrc_database = NULL; fm::bbrc_statistics = NULL; fm::bbrc_chisq = NULL; fm::bbrc_result = NULL;
      Reset();
      Defaults();
      SetType(_type);
      SetMinfreq(_minfreq);
      fm::bbrc_instance_present=true;
      if (getenv("FMINER_LAZAR")) fm::bbrc_do_yaml = false;
      if (getenv("FMINER_SMARTS")) fm::bbrc_gsp_out = false; 
      if (getenv("FMINER_PVALUES")) fm::bbrc_pvalues = true;
      if (getenv("FMINER_NO_AROMATIC_WC")) fm::bbrc_aromatic_wc = false;
      if (getenv("FMINER_SILENT")) {
        FILE* fp = freopen ("fminer_debug.txt","w",stderr);
      }
      if (getenv("FMINER_NR_HITS")) fm::bbrc_nr_hits = true;

  }
  else {
    cerr << "Error! Cannot create more than 1 instance." << endl; 
    exit(1);
  }

}

Bbrc::Bbrc(int _type, unsigned int _minfreq, float _chisq_val, bool _do_backbone) : init_mining_done(false) {
  if (!fm::bbrc_instance_present) {
      fm::bbrc_database = NULL; fm::bbrc_statistics = NULL; fm::bbrc_chisq = NULL; fm::bbrc_result = NULL;
      Reset();
      Defaults();
      SetType(_type);
      SetMinfreq(_minfreq);
      SetChisqSig(_chisq_val);
      SetBackbone(_do_backbone);
      fm::bbrc_instance_present=true;
      if (getenv("FMINER_LAZAR")) fm::bbrc_do_yaml = false;
      if (getenv("FMINER_SMARTS")) fm::bbrc_gsp_out = false; 
      if (getenv("FMINER_PVALUES")) fm::bbrc_pvalues = true;
      if (getenv("FMINER_NO_AROMATIC_WC")) fm::bbrc_aromatic_wc = false;
      if (getenv("FMINER_SILENT")) {
        FILE* fp = freopen ("fminer_debug.txt","w",stderr);
      }
      if (getenv("FMINER_NR_HITS")) fm::bbrc_nr_hits = true;

  }
  else {
    cerr << "Error! Cannot create more than 1 instance." << endl; 
    exit(1);
  }

}

Bbrc::~Bbrc() {
    if (fm::bbrc_instance_present) {
        delete fm::bbrc_database;
        delete fm::bbrc_statistics; 
        delete fm::bbrc_chisq; 
        delete fm::bbrc_ks;
        delete fm::bbrc_graphstate;
        delete fm::bbrc_closelegoccurrences;
        delete fm::bbrc_legoccurrences;

        fm::bbrc_Bbrccandidatelegsoccurrences.clear();
        fm::bbrc_candidatecloselegsoccs.clear();
        fm::bbrc_candidateBbrccloselegsoccsused.clear();

        fm::bbrc_instance_present=false;
    }
}

void Bbrc::Reset() { 
    if (fm::bbrc_instance_present) {
        delete fm::bbrc_database;
        delete fm::bbrc_statistics;
        delete fm::bbrc_chisq;
        delete fm::bbrc_ks;
        delete fm::bbrc_graphstate;
        delete fm::bbrc_closelegoccurrences;
        delete fm::bbrc_legoccurrences;
    }
    fm::bbrc_database = new BbrcDatabase();
    fm::bbrc_db_built = false;
    fm::bbrc_statistics = new BbrcStatistics();
    fm::bbrc_chisq = new ChisqBbrcConstraint(-1.0);
    fm::bbrc_ks = new KSBbrcConstraint(0.95);
    fm::bbrc_graphstate = new BbrcGraphState();
    fm::bbrc_closelegoccurrences = new CloseBbrcLegOccurrences();
    fm::bbrc_legoccurrences = new BbrcLegOccurrences();

    fm::bbrc_Bbrccandidatelegsoccurrences.clear();
    fm::bbrc_candidatecloselegsoccs.clear();
    fm::bbrc_candidateBbrccloselegsoccsused.clear();

    SetChisqActive(true); 
    fm::bbrc_result = &r;

    // clearing privates
    init_mining_done = false;
    comp_runner=0; 
    comp_no=0; 
    r.clear();
    inchi_compound_map.clear();
    inchi_compound_mmap.clear();
    activity_map.clear();
    weight_map.clear();

    if (getenv("FMINER_SILENT")) {
        fclose (stderr);
        FILE* fp = freopen ("fminer_debug.txt","w",stderr);
    }
}

void Bbrc::Defaults() {
    fm::bbrc_minfreq = 2;
    fm::bbrc_type = 2;
    fm::bbrc_do_pruning = true;
    fm::bbrc_console_out = false;
    fm::bbrc_aromatic = true;
    fm::bbrc_refine_singles = false;
    fm::bbrc_do_output=true;
    fm::bbrc_bbrc_sep=false;
    fm::bbrc_updated = true;
    fm::bbrc_gsp_out=true;
    fm::bbrc_nr_hits = false;

    // BBRC
    fm::bbrc_ks->sig = 0.95;
    fm::bbrc_chisq->sig = -1.0;
    fm::bbrc_do_backbone = true;
    fm::bbrc_adjust_ub = true;
    fm::bbrc_regression=false;
    fm::bbrc_do_yaml=true;
    fm::bbrc_pvalues=false;
    fm::bbrc_aromatic_wc=true;
}


// 2. Getter methods

int Bbrc::GetMinfreq(){return fm::bbrc_minfreq;}
int Bbrc::GetType(){return fm::bbrc_type;}
bool Bbrc::GetBackbone(){return fm::bbrc_do_backbone;}
bool Bbrc::GetDynamicUpperBound(){return fm::bbrc_adjust_ub;}
bool Bbrc::GetPruning() {return fm::bbrc_do_pruning;}
bool Bbrc::GetConsoleOut(){return fm::bbrc_console_out;}
bool Bbrc::GetAromatic() {return fm::bbrc_aromatic;}
bool Bbrc::GetRefineSingles() {return fm::bbrc_refine_singles;}
bool Bbrc::GetDoOutput() {return fm::bbrc_do_output;}
bool Bbrc::GetBbrcSep(){return fm::bbrc_bbrc_sep;}
bool Bbrc::GetChisqActive(){return fm::bbrc_chisq->active;}
float Bbrc::GetChisqSig(){if (!fm::bbrc_regression) return fm::bbrc_chisq->sig; else return fm::bbrc_ks->sig; }
bool Bbrc::GetRegression() {return fm::bbrc_regression;}



// 3. Setter methods

void Bbrc::SetMinfreq(int val) {
    // parameters not regarded in integrity constraints
    if (val < 1) { cerr << "Error! Invalid value '" << val << "' for parameter minfreq." << endl; exit(1); }
    if (val > 1 && GetRefineSingles()) { cerr << "Warning! Minimum frequency of '" << val << "' could not be set due to activated single refinement." << endl;}
    fm::bbrc_minfreq = val;
}

bool Bbrc::SetType(int val) {
    // parameters not regarded in integrity constraints
    if ((val != 1) && (val != 2)) { cerr << "Error! Invalid value '" << val << "' for parameter type." << endl; exit(1); }
    fm::bbrc_type = val;
    return 1;
}

bool Bbrc::SetBackbone(bool val) {
    // internal: chisq active
    if (val && !GetChisqActive()) {
        cerr << "Warning! BBRC mining could not be enabled due to deactivated significance criterium." << endl;
    }
    // -------- !db --------
    else if (!val && GetDynamicUpperBound()) {
        cerr << "Notice: Disabling dynamic upper bound pruning due to switched-off BBRC mining." << endl;
        SetDynamicUpperBound(false);
        fm::bbrc_do_backbone = val;
    }
    // -------- r!b ---------
    else if (val && GetBbrcSep()) {
        cerr << "Notice: Disabling BBRC separator due to enabled BBRC mining." << endl;
        SetBbrcSep(false);
        fm::bbrc_do_backbone = val;
    }
    else fm::bbrc_do_backbone = val;
    return 1;
}

bool Bbrc::SetDynamicUpperBound(bool val) {
    // -------- !db ---------
    if (val && !GetBackbone()) {
        cerr << "Warning! Dynamic upper bound pruning could not be enabled due to disabled BBRC mining." << endl;
    }
    // internal: chisq active
    else if (val && !GetChisqActive()) {
        cerr << "Warning! Dynamic upper bound pruning could not be enabled due to deactivated significance criterium." << endl;
    }
    // -------- !du ---------
    else if (val && !GetPruning()) {
        cerr << "Warning! Dynamic upper bound pruning could not be enabled due to deactivated statistical metric pruning." << endl;
    }
    else {
        fm::bbrc_adjust_ub=val; 
    }
    return 1;
}

bool Bbrc::SetPruning(bool val) {
    // internal: chisq active
    if (val && !GetChisqActive()) {
        cerr << "Warning! Statistical metric pruning could not be enabled due to deactivated significance criterium." << endl;
    }
    else {
        // ---------- !du -----------
        if (!val && GetDynamicUpperBound()) {
            cerr << "Notice: Disabling dynamic upper bound pruning due to disabled static upper bound pruning." << endl;
            SetDynamicUpperBound(false); 
        }
        // --------- ru -------
        if (!val && GetBbrcSep()) {
            cerr << "Notice: Disabling BBRC separator due to disabled static upper bound pruning." << endl;
            SetBbrcSep(false);
        }
        fm::bbrc_do_pruning=val;
    }
    return 1;
}

bool Bbrc::SetConsoleOut(bool val) {
    // console out not switched by fminer
    if (val) {
        if (GetBbrcSep()) cerr << "Warning! Console output could not be enabled due to enabled BBRC separator." << endl;
        else fm::bbrc_console_out=val;
    }
    return 1;
}

void Bbrc::SetAromatic(bool val) {
    fm::bbrc_aromatic = val;
}

bool Bbrc::SetRefineSingles(bool val) {
    fm::bbrc_refine_singles = val;
    // parameters not regarded in integrity constraints
    if (GetRefineSingles() && GetMinfreq() > 1) {
        cerr << "Notice: Using minimum frequency of 1 to refine singles." << endl;
        SetMinfreq(1);
    }
    return 1;
}

void Bbrc::SetDoOutput(bool val) {
    fm::bbrc_do_output = val;
}

bool Bbrc::SetBbrcSep(bool val) {
    //  ------- r!b ---------
    if (val && GetBackbone()) { 
        cerr << "Warning! BBRC separator could not be enabled due to enabled BBRC mining." << endl;
    }
    // -------- ru ---------
    if (val && !GetPruning()) {
        cerr << "Warning! BBRC separator could not be enabled due to disabled statistical metric pruning." << endl;
    }
    else {
        fm::bbrc_bbrc_sep=val;
        if (GetBbrcSep()) {
            // console out not switched by fminer
            if (GetConsoleOut()) {
                 cerr << "Notice: Disabling console output, using result vector." << endl;
                 SetConsoleOut(false);
            }
        }
    }
    return 1;
}

bool Bbrc::SetChisqActive(bool val) {
    fm::bbrc_chisq->active = val;
    // chisq active not switched by fminer
    if (!GetChisqActive()) {
        cerr << "Notice: Disabling dynamic upper bound pruning due to deactivated significance criterium." << endl;
        SetDynamicUpperBound(false); //order important
        cerr << "Notice: Disabling BBRC mining due to deactivated significance criterium." << endl;
        SetBackbone(false);
        cerr << "Notice: Disabling statistical metric pruning due to deactivated significance criterium." << endl;
        SetPruning(false);
        SetRegression(false);
    }
    return 1;
}

bool Bbrc::SetChisqSig(float _chisq_val) {
    // parameters not regarded in integrity constraints
    if (_chisq_val < 0.0 || _chisq_val > 1.0) { cerr << "Error! Invalid value '" << _chisq_val << "' for parameter chisq." << endl; exit(1); }
    if (fm::bbrc_regression) {
         fm::bbrc_ks->sig = _chisq_val; 
    }
    else {
        fm::bbrc_chisq->sig = gsl_cdf_chisq_Pinv(_chisq_val, 1);
        return 1;
    }
}

bool Bbrc::SetRegression(bool val) {
    fm::bbrc_regression = val;
    if (fm::bbrc_regression) {
         if (!GetBackbone()) {
            SetBackbone(true);
         }
         if (GetPruning()) {
            cerr << "Notice: Disabling statistical metric pruning due to activated regression." << endl;
            SetPruning(false);
         }
    }
    return 1;
}

// Forbidden in BBRC
bool Bbrc::SetMaxHops(int val) {
    return 0;
}

// 4. Other methods

vector<string>* Bbrc::MineRoot(unsigned int j) {
    fm::bbrc_result->clear();
    if (!init_mining_done) {
        if (!fm::bbrc_db_built) {
          AddDataCanonical();
        }
        // Adjust chisq bound
        if (!fm::bbrc_regression) {
          if (fm::bbrc_chisq->nr_acts.size()>1 && fm::bbrc_chisq->nr_acts.size() < 6) {
            if (fm::bbrc_chisq->sig == -1.0) { // do not override user-supplied threshold, only machine default.
              fm::bbrc_chisq->sig=fm::bbrc_chisq->df_thresholds[fm::bbrc_chisq->nr_acts.size()-1];
            }
          }
          else if (fm::bbrc_chisq->nr_acts.size()==1) {
            cout << "";
          }
          else {
            cerr << "Error! Too many classes: '" << fm::bbrc_chisq->nr_acts.size() << "' (Max. 5)." << endl;
            exit(1);
          }
        }
        fm::bbrc_database->edgecount (); 
        fm::bbrc_database->reorder (); 
        BbrcinitBbrcLegStatics (); 
        fm::bbrc_graphstate->init (); 
        if (fm::bbrc_bbrc_sep && !fm::bbrc_do_backbone && fm::bbrc_do_output && !fm::bbrc_console_out) (*fm::bbrc_result) << fm::bbrc_graphstate->sep();
        init_mining_done=true; 

        if (!fm::bbrc_regression) {
             cerr << "Settings:" << endl \
             << "---" << endl \
             << "Type:                                 " << GetType() << endl \
             << "Minimum frequency:                    " << GetMinfreq() << endl \
             << "Aromatic:                             " << GetAromatic() << endl \
             << "Chi-square active (chi-square-value): " << GetChisqActive() << " (" << GetChisqSig()<< ")" << endl \
             << "BBRC mining:                          " << GetBackbone() << endl \
             << "Statistical metric (dynamic) pruning: " << GetPruning() << " (" << GetDynamicUpperBound() << ")" << endl \
             << "Refine patterns with single support:  " << GetRefineSingles() << endl \
             << "Do output:                            " << GetDoOutput() << endl \
             << "BBRC sep:                             " << GetBbrcSep() << endl \
             << "Regression:                           " << GetRegression() << endl \
             << "---" << endl;
        }
        else {
             cerr << "Settings:" << endl \
             << "---" << endl \
             << "Type:                                 " << GetType() << endl \
             << "Minimum frequency:                    " << GetMinfreq() << endl \
             << "Aromatic:                             " << GetAromatic() << endl \
             << "KS active (p-value):                  " << GetChisqActive() << " (" << GetChisqSig()<< ")" << endl \
             << "BBRC mining:                          " << GetBackbone() << endl \
             << "Statistical metric (dynamic) pruning: " << GetPruning() << " (" << GetDynamicUpperBound() << ")" << endl \
             << "Refine patterns with single support:  " << GetRefineSingles() << endl \
             << "Do output:                            " << GetDoOutput() << endl \
             << "BBRC sep:                             " << GetBbrcSep() << endl \
             << "Regression:                           " << GetRegression() << endl \
             << "---" << endl;
        }

    }


    if (j >= fm::bbrc_database->nodelabels.size()) { cerr << "Error! Root node " << j << " does not exist." << endl;  exit(1); }
    if ( fm::bbrc_database->nodelabels[j].frequency >= fm::bbrc_minfreq && fm::bbrc_database->nodelabels[j].frequentedgelabels.size () ) {
        BbrcPath path(j);
        path.expand(); // mining step
    }
    if (getenv("FMINER_SILENT")) {
      fclose (stderr);
    }
    return fm::bbrc_result;
}

void Bbrc::ReadGsp(FILE* gsp){
    fm::bbrc_database->readGsp(gsp);
}

bool Bbrc::AddCompound(string smiles, unsigned int comp_id) {
  if (fm::bbrc_db_built) {
    cerr << "BbrcDatabase has been already processed! Please reset() and insert a new dataset." << endl;
    return false;
  }
  stringstream ss(smiles);
  OBConversion conv(&ss, &cout);
  if(!conv.SetInAndOutFormats("SMI","INCHI")) {
    cerr << "Formats not available" << endl;
    return false;
  }
  OBMol mol;
  if (!conv.Read(&mol)) {
    cerr << "Could not convert '" << smiles << "' (leaving out)." << endl;
    return false;
  }
  conv.SetOptions("w",OBConversion::OUTOPTIONS);
  string inchi = conv.WriteString(&mol);
  // remove newline
  string::size_type pos = inchi.find_last_not_of("\n");
  if (pos != string::npos) {
    inchi = inchi.substr(0, pos+1);
  }
  //cerr << "Inchi: '" << inchi << "'" << endl;

  // insert into map to check doubles
  pair<unsigned int, string> ori = make_pair(comp_id, smiles);
  pair< map<string,pair<unsigned int, string> >::iterator, bool> res = inchi_compound_map.insert(make_pair(inchi,ori));
  if (!res.second) {
    cerr << "Note: structure of '" << smiles << "' has been already inserted, inserting anyway..." << endl;
  }

  // insert into actual map augmented by number
  string inchi_no = inchi;
  inchi_no += "-";
  comp_runner++;
  stringstream out; out << comp_runner;
  string comp_runner_s = out.str();
  inchi_no += comp_runner_s;
  pair< map<string,pair<unsigned int, string> >::iterator, bool> resmm = inchi_compound_mmap.insert(make_pair(inchi_no,ori));
  return true;
}

bool Bbrc::AddActivity(float act, unsigned int comp_id) {
  if (fm::bbrc_db_built) {
    cerr << "BbrcDatabase has been already processed! Please reset() and insert a new dataset." << endl;
    return false;
  }
  activity_map.insert(make_pair(comp_id, act));
  return true;
}

bool Bbrc::AddWeight(float weight, unsigned int comp_id) {
  if (fm::bbrc_db_built) {
    cerr << "BbrcDatabase has been already processed! Please reset() and insert a new dataset." << endl;
    return false;
  }
  if (weight < 0.0) {
    cerr << "Weight '" << weight << "' for id '" << comp_id << "' is negative." << endl;
    return false;
  }
  weight_map.insert(make_pair(comp_id, weight));
  return true;
}



// the class factories
extern "C" Fminer* create0() {
    return new Bbrc();
}
extern "C" Fminer* create2(int _type, unsigned int _minfreq) {
    return new Bbrc(_type, _minfreq);
}
extern "C" Fminer* create4(int _type, unsigned int _minfreq, float _chisq_val, bool _do_backbone) {
    return new Bbrc(_type, _minfreq, _chisq_val, _do_backbone);
}
extern "C" void destroy(Fminer* f) {
    delete f;
}

extern "C" void usage() {
    cerr << endl;
    cerr << "Options for Usage 1 (BBRC mining using dynamic upper bound pruning): " << endl;
    cerr << "       [-f minfreq] [-l type] [-s] [-a] [-o] [-g] [-d [-b [-u]]] [-p p_value]" << endl;
    cerr << endl;
    cerr << "Options for Usage 2 (Frequent subgraph mining): " << endl;
    cerr << "       [-f minfreq] [-l type] [-s] [-a] [-o]" << endl;
    cerr << endl;
}

bool Bbrc::AddDataCanonical() {
    // AM: now insert all structures into the database
    // in canonical ordering according to inchis
    comp_runner=0;
    if (fm::bbrc_regression) {
      vector<float> activity_values;
      for (map<string, pair<unsigned int, string> >::iterator it = inchi_compound_mmap.begin(); it != inchi_compound_mmap.end(); it++) {
        if (activity_map.find(it->second.first) != activity_map.end()) { // need this for reading value in the next line!
          activity_values.push_back(activity_map.find(it->second.first)->second); // 0 gets push w/o previous check line!
        }
      }
      // cut quantiles
      float min_thr, max_thr = 0.0;
      min_thr = quantile(activity_values,0.0125); // find 1.25% quantile
      max_thr = quantile(activity_values,0.9875); // find 98.75% quantile

      for (map<string, pair<unsigned int, string> >::iterator it = inchi_compound_mmap.begin(); it != inchi_compound_mmap.end(); it++) {
        if (activity_map.find(it->second.first) != activity_map.end()) {
          activity_map[it->second.first]=activity_map[it->second.first];
          if (activity_map[it->second.first] < min_thr) activity_map[it->second.first]=min_thr;
          if (activity_map[it->second.first] > max_thr) activity_map[it->second.first]=max_thr;
        }
      }
    }

    // AM: weight map initialization to 1 for all instances
    if (weight_map.size() == 0) { // when user has done nothing
      for (map<string, pair<unsigned int, string> >::iterator it = inchi_compound_mmap.begin(); it != inchi_compound_mmap.end(); it++) {
        weight_map.insert(make_pair(it->second.first, 1.0));
      }
    }

    for (map<string, pair<unsigned int, string> >::iterator it = inchi_compound_mmap.begin(); it != inchi_compound_mmap.end(); it++) {
      AddCompoundCanonical(it->second.second, it->second.first); // smiles, comp_id
      float activity = activity_map.find(it->second.first)->second;
      AddActivityCanonical(activity, it->second.first); // act, comp_id
    }

    for (map<unsigned int, float>::iterator it = weight_map.begin(); it != weight_map.end(); it++) {
      float weight = it->second;
      unsigned int comp_id = it->first;
      CheckWeight(weight, comp_id); // weight, comp_id: remove weights of non-existing structures from map
    }
    NormalizeWeights(weight_map);

    for (map<string, pair<unsigned int, string> >::iterator it = inchi_compound_mmap.begin(); it != inchi_compound_mmap.end(); it++) {
      float weight = weight_map.find(it->second.first)->second;
      AddWeightCanonical(weight, it->second.first); // weight, comp_id
    }

    fm::bbrc_db_built=true;
    inchi_compound_map.clear();
    inchi_compound_mmap.clear();
    activity_map.clear();
    weight_map.clear();
}

bool Bbrc::AddCompoundCanonical(string smiles, unsigned int comp_id) {
  bool insert_done=false;
  if (comp_id<=0) { cerr << "Error! IDs must be of type: Int > 0." << endl;}
  else {
    if (activity_map.find(comp_id) == activity_map.end() && GetChisqActive()) {
      cerr << "Error on compound '" << comp_runner << "', id '" << comp_id << "': no activity found." << endl;
      return false;
    }
    if (weight_map.find(comp_id) == weight_map.end() && GetChisqActive()) {
      cerr << "Error on compound '" << comp_runner << "', id '" << comp_id << "': no weight found." << endl;
      return false;
    }
    else {
      if (fm::bbrc_database->readTreeSmi (smiles, comp_no, comp_id, comp_runner)) {
        insert_done=true;
        comp_no++;
      }
      else { cerr << "Error on compound '" << comp_runner << "', id '" << comp_id << "'." << endl; }
    }
    comp_runner++;
  }
  return insert_done;
}

bool Bbrc::AddActivityCanonical(float act, unsigned int comp_id) {
  if (fm::bbrc_database->trees_map.find(comp_id) == fm::bbrc_database->trees_map.end()) { 
    cerr << "No structure for ID " << comp_id << " when adding activity. Ignoring entry!" << endl; return false; 
  }
  else {
    if (!fm::bbrc_regression) {
      AddChiSq(fm::bbrc_database->trees_map[comp_id]->activity = act);
    }
    else {
      AddKS(fm::bbrc_database->trees_map[comp_id]->activity = act);
    }
    return true;
  }
}

bool Bbrc::CheckWeight(float weight, unsigned int comp_id) {
  if (fm::bbrc_database->trees_map.find(comp_id) == fm::bbrc_database->trees_map.end()) { 
    weight_map.erase(comp_id);
    cout << "FPP" << endl;
    return false;
  }
  return true;
}

bool Bbrc::NormalizeWeights(map<unsigned int, float> weight_map) {
  map<unsigned int, float>::iterator weight_map_it;
  float weight_sum = 0.0;
  for (weight_map_it = weight_map.begin(); weight_map_it != weight_map.end(); weight_map_it++) {
    cout << "AM: w " << weight_map_it->second << endl;
  }
  for (weight_map_it = weight_map.begin(); weight_map_it != weight_map.end(); weight_map_it++) {
    weight_sum += weight_map_it->second;
  }
  int nr_weights = weight_map.size();
  for (weight_map_it = weight_map.begin(); weight_map_it != weight_map.end(); weight_map_it++) {
    weight_map_it->second = weight_map_it->second * nr_weights;
    weight_map_it->second = weight_map_it->second / weight_sum;
  }
  cout << endl;
  for (weight_map_it = weight_map.begin(); weight_map_it != weight_map.end(); weight_map_it++) {
    cout << "AM: w " << weight_map_it->second << endl;
  }
  cout << endl;
}

bool Bbrc::AddWeightCanonical(float weight, unsigned int comp_id) {
  if (fm::bbrc_database->trees_map.find(comp_id) == fm::bbrc_database->trees_map.end()) { 
    cerr << "No structure for ID " << comp_id << " when adding weight. Ignoring entry!" << endl; return false; 
  }
  else {
    if (!fm::bbrc_regression) {
      fm::bbrc_database->trees_map[comp_id]->weight = weight;
    }
    return true;
  }
}
