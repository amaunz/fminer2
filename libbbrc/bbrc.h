// bbrc.h
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

#ifndef BBRC_H
#define BBRC_H

#include "../fminer/fminer.h"
#include "misc.h"
#include "closeleg.h"
#include "graphstate.h"

namespace fm { 

    extern bool bbrc_adjust_ub;
    extern bool bbrc_do_pruning;
    extern bool bbrc_aromatic;
    extern ChisqBbrcConstraint* bbrc_chisq;
    extern KSBbrcConstraint* bbrc_ks;
    extern bool bbrc_do_yaml;
    extern bool bbrc_gsp_out;
    extern bool bbrc_bbrc_sep;
    extern bool bbrc_regression;
    extern bool bbrc_db_built;

}

class Bbrc : public Fminer {

  public:

    /** @name Inits
     *  Initializer functions.
     */
    //@{
    Bbrc (); //!< Constructor for standard settings: 95% significance Bbrclevel, minimum frequency 2, type trees, dynamic upper bound, BBRC.
    Bbrc (int _type, unsigned int _minfreq); //!< Like standard constructor, but type and minimum frequency configurable.
    Bbrc (int _type, unsigned int _minfreq, float chisq_val, bool _do_backbone); //!< Like standard constructor, but type, minimum frequency, significance Bbrclevel and BBRC configurable.

    ~Bbrc();
    void Reset(); //!< Use this to clear the database before feeding new compounds and activities.
    void Defaults(); //!< Use this to set default parameters as in default constructor.
    //@}

    /** @name Getters
     *  Getter functions.
     */
    //@{
    int GetMinfreq(); //!< Get minimum frequency.
    int GetType(); //!< Get type.
    bool GetBackbone(); //!< Get whether BBRC representatives should be mined.
    bool GetDynamicUpperBound(); //!< Get whether dynamic upper bound pruning is used.
    bool GetPruning(); //!< Get whether statistical metric pruning should be used.
    bool GetConsoleOut(); //!< Get whether output should be directed to the console.
    bool GetAromatic(); //!< Get whether aromatic rings should be perceived instead of Kekule notation.
    bool GetRefineSingles(); //!< Get whether fragments with frequency 1 should be refined.
    bool GetDoOutput(); //!< Get whether output is enabled.
    bool GetBbrcSep(); //!< Get whether BBRCs should be separated in the output.
    bool GetChisqActive(); //!< Get whether chi-square filter is active.
    float GetChisqSig(); //!< Get significance threshold.
    bool GetRegression(); //!< Get whether continuous activity values should be used.
    //@}

    /** @name Setters
     *  Setter functions.
     */
    //@{
    void SetMinfreq(int val); //!< Set minimum frequency (>=1 here). Same as '-f'.
    bool SetType(int val); //!< Set type 1 (paths) or 2 (trees) here. Same as '-l'.
    bool SetBackbone(bool val); //!< Pass 'false' here to switch off mining for BBRC representatives. Same as '-b'.
    bool SetDynamicUpperBound(bool val); //!< Pass 'false' here to disable dynamic upper bound pruning (e.g. for performance measures). Same as '-d'.
    bool SetPruning(bool val); //!< Pass 'false' here to disable statistical metrical pruning completely. Same as '-u'.
    bool SetConsoleOut(bool val); //!< Pass 'true' here to disable usage of result vector and directly print each fragment to the console (saves memory).
    void SetAromatic(bool val); //!< Pass 'true' here to enable aromatic rings and use Kekule notation. IMPORTANT! SET THIS BEFORE CALLING AddCompound()! Same as '-a'.
    bool SetRefineSingles(bool val); //!< Pass 'true' here to enable refinement of fragments with frequency 1. Same as '-s'.
    void SetDoOutput(bool val); //!< Pass 'false' here to disable output. Same as '-o'.
    bool SetBbrcSep(bool val); //!< Set this to 'true' to enable BBRC separators in output.
    bool SetChisqActive(bool val); //!< Set this to 'true' to enable chi-square filter.
    bool SetChisqSig(float _chisq_val); //!< Set significance threshold here (between 0 and 1). Same as '-p'.
    bool SetRegression(bool val); //!< Set 'true' here to enable continuous activity values. Same as '-g'.
    bool SetMaxHops(int val); //!< Dummy method for max hops (only used in LAST-PM). Same as '-m'.
    //@}
    /** @name Others
     *  Other functions.
     */
    //@{
    vector<string>* MineRoot(unsigned int j); //!< Mine fragments rooted at the j-th root node (element type).
    void ReadGsp(FILE* gsp); //!< Read in a gSpan file
    bool AddCompound(string smiles, unsigned int comp_id); //!< Add a compound to the database.
    // KS: bool AddActivity(bool act, unsigned int comp_id); //!< Add an activity to the database.
    // KS: recognize regr field
    bool AddActivity(float act, unsigned int comp_id); //!< Add an activity to the database.
    bool AddWeight(float weight, unsigned int comp_id); //!< Add a weight to the database.
    int GetNoRootNodes() {if (!fm::bbrc_db_built) AddDataCanonical() ; return fm::bbrc_database->nodelabels.size();} //!< Get number of root nodes (different element types).
    int GetNoCompounds() {if (!fm::bbrc_db_built) AddDataCanonical() ; return fm::bbrc_database->trees.size();} //!< Get number of compounds in the database.
    //@}
    
  private:
    void AddChiSq(float act){fm::bbrc_chisq->nr_acts[act]++; fm::bbrc_chisq->n++;}
    // KS: Insert value into set of activities
    void AddKS(float val){fm::bbrc_ks->all.push_back(val);}

    bool init_mining_done;
    int comp_runner;
    int comp_no;

    vector<string> r;
    // ONLY FOR INTERNAL USE. DO NOT MAKE PUBLIC!
    map<string, pair<unsigned int, string> > inchi_compound_map;    // AM: structure inchi => (id, smi) for canonical input to check for double structures
    map<string, pair<unsigned int, string> > inchi_compound_mmap;   // AM: structure inchi => (id, smi) for canonical input to use for actual storage
    map<unsigned int, float> activity_map;                          // AM: structure id => act for canonical input
    map<unsigned int, float> weight_map;                            // AM: structure id => weight for canonical input

    bool AddDataCanonical();                                        //!< Only to be called by MineRoot!
    bool AddCompoundCanonical(string smiles, unsigned int comp_id); //!< Only to be called by AddDataCanonical!
    bool AddActivityCanonical(float act, unsigned int comp_id);     //!< Only to be called by AddDataCanonical!
    bool CheckWeight(float weight, unsigned int comp_id);           //!< Only to be called by AddDataCanonical!
    bool NormalizeWeights(map<unsigned int, float> weight_map);     //!< Only to be called by AddDataCanonical!
    bool AddWeightCanonical(float weight, unsigned int comp_id);    //!< Only to be called by AddDataCanonical!

};

#endif
