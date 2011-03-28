// last.cpp
// modified 2010
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

#include "last.h"
#include "globals.h"


// 1. Constructors and Initializers

Last::Last() : init_mining_done(false) {
  if (!fm::last_instance_present) {
      fm::last_database = NULL; fm::last_statistics = NULL; fm::last_chisq = NULL; fm::last_result = NULL;
      Reset();
      Defaults();
      fm::last_instance_present=true;
      fm::last_gsp_out = false; 
      if (getenv("FMINER_SILENT")) {
        FILE* fp = freopen ("fminer_last_debug.txt","w",stderr);
      }
  }
  else {
    cerr << "Error! Cannot create more than 1 instance." << endl; 
    exit(1);
  }
}

/*
Last::Last(int _type, unsigned int _minfreq) : init_mining_done(false) {
  if (!fm::last_instance_present) {
      fm::last_database = NULL; fm::last_statistics = NULL; fm::last_chisq = NULL; fm::last_result = NULL;
      Reset();
      Defaults();
      SetType(_type);
      SetMinfreq(_minfreq);
      fm::last_instance_present=true;
      fm::last_gsp_out = false; 
  }
  else {
    cerr << "Error! Cannot create more than 1 instance." << endl; 
    exit(1);
  }
}

Last::Last(int _type, unsigned int _minfreq, float _chisq_val, bool _do_backbone) : init_mining_done(false) {
  if (!fm::last_instance_present) {
      fm::last_database = NULL; fm::last_statistics = NULL; fm::last_chisq = NULL; fm::last_result = NULL;
      Reset();
      Defaults();
      SetType(_type);
      SetMinfreq(_minfreq);
      SetChisqSig(_chisq_val);
      fm::last_instance_present=true;
      fm::last_gsp_out = false; 

  }
  else {
    cerr << "Error! Cannot create more than 1 instance." << endl; 
    exit(1);
  }
}
*/

Last::~Last() {
    if (fm::last_instance_present) {
        delete fm::last_database;
        delete fm::last_statistics; 
        delete fm::last_chisq; 
        delete fm::last_ks;
        delete fm::last_graphstate;
        delete fm::last_closelegoccurrences;
        delete fm::last_legoccurrences;

        fm::last_Lastcandidatelegsoccurrences.clear();
        fm::last_candidatecloselegsoccs.clear();
        fm::last_candidateLastcloselegsoccsused.clear();

        fm::last_instance_present=false;
    }
}

void Last::Reset() { 
    if (fm::last_instance_present) {
        delete fm::last_database;
        delete fm::last_statistics;
        delete fm::last_chisq;
        delete fm::last_ks;
        delete fm::last_graphstate;
        delete fm::last_closelegoccurrences;
        delete fm::last_legoccurrences;
    }
    fm::last_database = new LastDatabase();
    fm::last_db_built = false;
    fm::last_statistics = new LastStatistics();
    fm::last_chisq = new ChisqLastConstraint(3.84146);
    fm::last_ks = new KSLastConstraint(0.95);
    fm::last_graphstate = new LastGraphState();
    fm::last_closelegoccurrences = new CloseLastLegOccurrences();
    fm::last_legoccurrences = new LastLegOccurrences();

    fm::last_candidateLastcloselegsoccsused.clear();

    fm::last_chisq->active=true; 
    fm::last_result = &r;
    comp_runner=0; 
    comp_no=0; 
    init_mining_done = false;
}

void Last::Defaults() {
    fm::last_minfreq = 2;
    fm::last_type = 2;
    fm::last_do_pruning = true;
    fm::last_console_out = true;
    fm::last_aromatic = true;
    fm::last_refine_singles = false;
    fm::last_do_output=true;
    fm::last_bbrc_sep=false;
    fm::last_updated = true;
    fm::last_gsp_out=true;
    fm::last_regression=false;

    // LAST
    fm::last_do_last=true;
    fm::last_hops=0;
    fm::last_die = 0;
    fm::last_max_hops = 25;
}


// 2. Getter methods

int Last::GetMinfreq(){return fm::last_minfreq;}
int Last::GetType(){return fm::last_type;}
bool Last::GetBackbone(){return false;}
bool Last::GetDynamicUpperBound(){return false;}
bool Last::GetPruning() {return fm::last_do_pruning;}
bool Last::GetConsoleOut(){return fm::last_console_out;}
bool Last::GetAromatic() {return fm::last_aromatic;}
bool Last::GetRefineSingles() {return fm::last_refine_singles;}
bool Last::GetDoOutput() {return fm::last_do_output;}
bool Last::GetBbrcSep(){return fm::last_bbrc_sep;}
bool Last::GetChisqActive(){return fm::last_chisq->active;}
float Last::GetChisqSig(){if (!fm::last_regression) return fm::last_chisq->sig; else return fm::last_ks->sig;}
bool Last::GetRegression() {return fm::last_regression;}
int Last::GetMaxHops() {return fm::last_max_hops;}


// 3. Setter methods

void Last::SetMinfreq(int val) {
    if (val < 1) { cerr << "Error! Invalid value '" << val << "' for parameter minfreq." << endl; exit(1); }
    if (val > 1 && GetRefineSingles()) { cerr << "Warning! Minimum frequency of '" << val << "' could not be set due to activated single refinement." << endl;}
    fm::last_minfreq = val;
}

// These methods report forbidden switches (synopsis) back to main
// They also report forbidden argument switches (exception: arguments equal defaults)
bool Last::SetType(int val) {
    return 0;
}

bool Last::SetBackbone(bool val) {
    return 0;
}

bool Last::SetDynamicUpperBound(bool val) {
    return 0;
}

bool Last::SetPruning(bool val) {
    return 0;
}

bool Last::SetConsoleOut(bool val) {
    // console out not switched by fminer
    fm::last_console_out=val;
    return 1;
}

void Last::SetAromatic(bool val) {
    fm::last_aromatic = val;
}

bool Last::SetRefineSingles(bool val) {
    return 0;
}

void Last::SetDoOutput(bool val) {
    fm::last_do_output = val;
}

bool Last::SetBbrcSep(bool val) {
    return 0;
}

bool Last::SetChisqActive(bool val) {
    return 0;
}

bool Last::SetChisqSig(float _chisq_val) {
    return 0;
}

bool Last::SetRegression(bool val) {
    // return 0;
    // TODO: enable regression
    fm::last_regression=val;
    return 1;
}

bool Last::SetMaxHops(int val) {
    fm::last_max_hops=val;
    return 1;
}

// 4. Other methods

vector<string>* Last::MineRoot(unsigned int j) {
    fm::last_result->clear();
    if (!init_mining_done) {
        if (!fm::last_db_built) {
          AddDataCanonical();
        }
        if (fm::last_chisq->active) {
            each (fm::last_database->trees) {
                if (fm::last_database->trees[i]->activity == -1) {
                    cerr << "Error! ID " << fm::last_database->trees[i]->orig_tid << " is missing activity information." << endl;
                    exit(1);
                }
            }
        }
        fm::last_database->edgecount (); 
        fm::last_database->reorder (); 
        LastinitLastLegStatics (); 
        fm::last_graphstate->init (); 
        if (fm::last_bbrc_sep && fm::last_do_output && !fm::last_console_out) (*fm::last_result) << fm::last_graphstate->sep();
        init_mining_done=true; 

        if (!fm::last_regression) {
            cerr << "Settings:" << endl \
                 << "---" << endl \
                 << "Type:                                 " << GetType() << endl \
                 << "Minimum frequency:                    " << GetMinfreq() << endl \
                 << "Aromatic:                             " << GetAromatic() << endl \
                 << "Regression:                           " << GetRegression() << endl \
                 << "Chi-square active (chi-square-value): " << GetChisqActive() << " (" << GetChisqSig()<< ")" << endl \
                 << "Statistical metric pruning:           " << GetPruning() << endl \
                 << "Do output:                            " << GetDoOutput() << endl \
                 << "Max Hops:                             " << GetMaxHops() << endl \
                 << "---" << endl;
        }
        else {
            cerr << "Settings:" << endl \
                 << "---" << endl \
                 << "Type:                                 " << GetType() << endl \
                 << "Minimum frequency:                    " << GetMinfreq() << endl \
                 << "Aromatic:                             " << GetAromatic() << endl \
                 << "Regression:                           " << GetRegression() << endl \
                 << "KS active (p-value):                  " << GetChisqActive() << " (" << GetChisqSig()<< ")" << endl \
                 << "Statistical metric pruning:           " << GetPruning() << endl \
                 << "Do output:                            " << GetDoOutput() << endl \
                 << "Max Hops:                             " << GetMaxHops() << endl \
                 << "---" << endl;
        }



        if (fm::last_do_output) {

          string xml_header = 
"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n\
xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n\
xsi:noNamespaceSchemaLocation=\"graphml.xsd\">\n\
\n\
<!-- LAtent STructure Mining (LAST) descriptors-->\n\
\n\
<key id=\"act\" for=\"graph\" attr.name=\"activating\" attr.type=\"boolean\" />\n\
<key id=\"hops\" for=\"graph\" attr.name=\"hops\" attr.type=\"int\" />\n\
<key id=\"lab_n\" for=\"node\" attr.name=\"node_labels\" attr.type=\"string\" />\n\
<key id=\"lab_e\" for=\"edge\" attr.name=\"edge_labels\" attr.type=\"string\" />\n\
<key id=\"weight\" for=\"edge\" attr.name=\"edge_weight\" attr.type=\"int\" />\n\
<key id=\"del\" for=\"edge\" attr.name=\"edge_deleted\" attr.type=\"boolean\" />\n\n";

         if (!fm::last_console_out) (*fm::last_result) << xml_header;
         else cout << xml_header;
      }
    }

    if (j >= fm::last_database->nodelabels.size()) { cerr << "Error! Root node " << j << " does not exist." << endl;  exit(1); }
    if ( fm::last_database->nodelabels[j].frequency >= fm::last_minfreq && fm::last_database->nodelabels[j].frequentedgelabels.size () ) {
        LastPath path(j);
        path.expand(); // mining step
    }
    if (j==GetNoRootNodes()-1 && fm::last_do_output) {
      if (!fm::last_console_out) (*fm::last_result) << "</graphml>\n";
      else cout << "</graphml>" << endl;
    }
    if (getenv("FMINER_SILENT")) {
      fclose (stderr);
    }
    return fm::last_result;
}

void Last::ReadGsp(FILE* gsp){
    fm::last_database->readGsp(gsp);
}

bool Last::AddCompound(string smiles, unsigned int comp_id) {
  if (fm::last_db_built) {
    cerr << "LastDatabase has been already processed! Please reset() and insert a new dataset." << endl;
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
  }
  conv.SetOptions("w",OBConversion::OUTOPTIONS);
  string inchi = conv.WriteString(&mol);
  // remove newline
  string::size_type pos = inchi.find_last_not_of("\n");
  if (pos != string::npos) {
    inchi = inchi.substr(0, pos+1);
  }
  //cerr << "Inchi: '" << inchi << "'" << endl;
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

bool Last::AddActivity(float act, unsigned int comp_id) {
  if (fm::last_db_built) {
    cerr << "LastDatabase has been already processed! Please reset() and insert a new dataset." << endl;
    return false;
  }
  activity_map.insert(make_pair(comp_id, act));
  return true;
}



// the class factories
extern "C" Fminer* create0() {
    return new Last();
}
/*
extern "C" Fminer* create2(int _type, unsigned int _minfreq) {
    return new Last(_type, _minfreq);
}
extern "C" Fminer* create4(int _type, unsigned int _minfreq, float _chisq_val, bool _do_backbone) {
    return new Last(_type, _minfreq, _chisq_val, _do_backbone);
}
*/
extern "C" void destroy(Fminer* l) {
    delete l;
}

extern "C" void usage() {
    cerr << endl;
    cerr << "Options for Usage 1 (LAtent STructure-Pattern Mining): " << endl;
    cerr << "       [-f minfreq] [-m maxhops] [-a] [-o] [-g] <graphs> <activities>" << endl;
    cerr << endl;
}

bool Last::AddDataCanonical() {
    // AM: now insert all structures into the database
    // in canonical ordering according to inchis
    comp_runner=0;
    for (map<string, pair<unsigned int, string> >::iterator it = inchi_compound_mmap.begin(); it != inchi_compound_mmap.end(); it++) {
      //cerr << it->second.first << "\t" << it->second.second << endl;
      AddCompoundCanonical(it->second.second, it->second.first); // smiles, comp_id
      AddActivityCanonical(activity_map[it->second.first], it->second.first); // act, comp_id
    }
    fm::last_db_built=true;
    inchi_compound_map.clear();
    inchi_compound_mmap.clear();
}

bool Last::AddCompoundCanonical(string smiles, unsigned int comp_id) {
  bool insert_done=false;
  if (comp_id<=0) { cerr << "Error! IDs must be of type: Int > 0." << endl;}
  else {
    if (fm::last_database->readTreeSmi (smiles, comp_no, comp_id, comp_runner)) {
      insert_done=true;
      comp_no++;
    }
    else { cerr << "Error on compound " << comp_runner << ", id " << comp_id << "." << endl; }
    comp_runner++;
  }
  return insert_done;
}

bool Last::AddActivityCanonical(float act, unsigned int comp_id) {
  if (fm::last_database->trees_map[comp_id] == NULL) { 
    cerr << "No structure for ID " << comp_id << ". Ignoring entry!" << endl; return false; 
  }
  else {
    if (!fm::last_regression) {
      if ((fm::last_database->trees_map[comp_id]->activity = act) == 1.0) AddChiSqNa();
      else AddChiSqNi();
    }
    else {
      if ((fm::last_database->trees_map[comp_id]->activity = act)) AddKS(act);
    }
    return true;
  }
}

float Last::ChisqTest(vector<float> all, vector<float> feat) {
  int x=0; int y=0; int m=0; int n=0;
  each(all) {
    if ((all[i] != 0.0) && (all[i] != 1.0)) {
      cerr << "liblast: got invalid values for chi-square test." << endl;
      exit(1);
    }
    if (all[i] == 1.0) m++;
    n++;
  }
  each(feat) {
    if ((feat[i] != 0.0) && (feat[i] != 1.0)) {
      cerr << "liblast: got invalid values for chi-square test." << endl;
      exit(1);
    }
    if (feat[i] == 1.0) y++;
    x++;
  }
  return fm::last_chisq->ChiSqTest(x,y,m,n-m);
}
