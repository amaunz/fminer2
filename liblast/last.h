// last.h
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

#ifndef LAST_H
#define LAST_H

#include "../fminer/fminer.h"
#include "misc.h"
#include "closeleg.h"
#include "graphstate.h"

namespace fm { 

    extern bool do_pruning;
    extern bool aromatic;
    extern ChisqConstraint* chisq;
    extern KSConstraint* ks;
    extern bool gsp_out;
    extern bool bbrc_sep;
    extern bool regression;
    extern int max_hops;
    extern bool db_built;

}

class Last : public Fminer {

  public:

    /** @name Inits
     *  Initializer functions.
     */
    //@{
    Last (); //!< Constructor for standard settings: 95% significance level, minimum frequency 2, type trees, dynamic upper bound, BBRC.
    /*
    Last (int _type, unsigned int _minfreq); //!< Like standard constructor, but type and minimum frequency configurable.
    Last (int _type, unsigned int _minfreq, float chisq_val, bool _do_backbone); //!< Like standard constructor, but type, minimum frequency, significance level and BBRC configurable.
    */
    ~Last();
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
    bool GetRegression(); //!< Dummy method for regression (only used for bbrcs).
    int GetMaxHops(); //!< Get maximum number of hops.

    //@}

    /** @name Setters
     *  Setter functions.
     */
    //@{
    void SetMinfreq(int val); //!< Set minimum frequency (>=1 here).
    bool SetType(int val); //!< Set type 1 (paths) or 2 (trees) here.
    bool SetBackbone(bool val); //!< Pass 'false' here to switch off mining for BBRC representatives.
    bool SetDynamicUpperBound(bool val); //!< Pass 'false' here to disable dynamic upper bound pruning (e.g. for performance measures).
    bool SetPruning(bool val); //!< Pass 'false' here to disable statistical metrical pruning completely.
    bool SetConsoleOut(bool val); //!< Pass 'true' here to disable usage of result vector and directly print each fragment to the console (saves memory).
    void SetAromatic(bool val); //!< Pass 'true' here to enable aromatic rings and use Kekule notation.
    bool SetRefineSingles(bool val); //!< Pass 'true' here to enable refinement of fragments with frequency 1.
    void SetDoOutput(bool val); //!< Pass 'false' here to disable output.
    bool SetBbrcSep(bool val); //!< Set this to 'true' to enable BBRC separators in output.
    bool SetChisqActive(bool val); //!< Set this to 'true' to enable chi-square filter.
    bool SetChisqSig(float _chisq_val); //!< Set significance threshold here (between 0 and 1).
    bool SetRegression(bool val); //!< Dummy method for regression (only used for bbrcs).
    bool SetMaxHops(int val); //!< Set maximum number of hops.
    //@}
    
    /** @name Others
     *  Other functions.
     */
    //@{
    vector<string>* MineRoot(unsigned int j); //!< Mine fragments rooted at the j-th root node (element type).
    void ReadGsp(FILE* gsp); //!< Read in a gSpan file
    bool AddCompound(string smiles, unsigned int comp_id); //!< Add a compound to the database.
    bool AddActivity(float act, unsigned int comp_id); //!< Add an activity to the database.
    int GetNoRootNodes() {if (!fm::db_built) AddDataCanonical() ; return fm::database->nodelabels.size();} //!< Get number of root nodes (different element types).
    int GetNoCompounds() {if (!fm::db_built) AddDataCanonical() ; return fm::database->trees.size();} //!< Get number of compounds in the database.
    //@}
    
  private:
    void AddChiSqNa(){fm::chisq->na++;fm::chisq->n++;}
    void AddChiSqNi(){fm::chisq->ni++;fm::chisq->n++;}
    // KS: Insert value into set of activities
    void AddKS(float val){fm::ks->all.push_back(val);}

    bool init_mining_done;
    int comp_runner;
    int comp_no;

    vector<string> r;
    // ONLY FOR INTERNAL USE. DO NOT MAKE PUBLIC!
    map<string, pair<unsigned int, string> > inchi_compound_map; // AM: structure inchi => (id, smi) for canonical input
    map<unsigned int, float> activity_map;                       // AM: structure inchi => (id, smi) for canonical input
    bool AddDataCanonical();                                        //!< Only to be called by MineRoot!
    bool AddCompoundCanonical(string smiles, unsigned int comp_id); //!< Only to be called by AddDataCanonical!
    bool AddActivityCanonical(float act, unsigned int comp_id);     //!< Only to be called by AddDataCanonical!

};

#endif
