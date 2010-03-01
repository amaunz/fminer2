// last.h
// modified 2010
// © 2008 by Andreas Maunz, andreas@maunz.de, oct 2008

/*
    This file is part of LibFminer (libfminer).

    LibFminer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LibFminer is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LibFminer.  If not, see <http://www.gnu.org/licenses/>.
 */

//#ifndef LAST_H
//#define LAST_H

#include "../fminer/fminer.h"
#include "misc.h"
#include "closeleg.h"
#include "graphstate.h"

namespace fm { 

    extern bool do_pruning;
    extern bool aromatic;
    extern ChisqConstraint* chisq;
    extern bool gsp_out;
    extern bool bbrc_sep;
    extern bool line_nrs;

}
class Last : public Fminer {

  public:

    /** @name Inits
     *  Initializer functions.
     */
    //@{
    Last (); //!< Constructor for standard settings: 95% significance level, minimum frequency 2, type trees, dynamic upper bound, BBRC.
    Last (int _type, unsigned int _minfreq); //!< Like standard constructor, but type and minimum frequency configurable.
    Last (int _type, unsigned int _minfreq, float chisq_val, bool _do_backbone); //!< Like standard constructor, but type, minimum frequency, significance level and BBRC configurable.

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
    bool GetMostSpecTreesOnly(); //!< Get whether most specific trees only should be mined for every BBRC.
    bool GetChisqActive(); //!< Get whether chi-square filter is active.
    float GetChisqSig(); //!< Get significance threshold.
    bool GetLineNrs(); //!< Get whether line numbers should be used in the output file.
    bool GetRegression(); //!< Dummy method for regression (only used for bbrcs).

    //@}

    /** @name Setters
     *  Setter functions.
     */
    //@{
    void SetMinfreq(int val); //!< Set minimum frequency (>=1 here).
    void SetType(int val); //!< Set type 1 (paths) or 2 (trees) here.
    void SetBackbone(bool val); //!< Pass 'false' here to switch off mining for BBRC representatives.
    void SetDynamicUpperBound(bool val); //!< Pass 'false' here to disable dynamic upper bound pruning (e.g. for performance measures).
    void SetPruning(bool val); //!< Pass 'false' here to disable statistical metrical pruning completely.
    void SetConsoleOut(bool val); //!< Pass 'true' here to disable usage of result vector and directly print each fragment to the console (saves memory).
    void SetAromatic(bool val); //!< Pass 'true' here to enable aromatic rings and use Kekule notation.
    void SetRefineSingles(bool val); //!< Pass 'true' here to enable refinement of fragments with frequency 1.
    void SetDoOutput(bool val); //!< Pass 'false' here to disable output.
    void SetBbrcSep(bool val); //!< Set this to 'true' to enable BBRC separators in output.
    void SetMostSpecTreesOnly(bool val); //!< Set this to 'true' to enable mining for the most specific tree patterns only.
    void SetChisqActive(bool val); //!< Set this to 'true' to enable chi-square filter.
    void SetChisqSig(float _chisq_val); //!< Set significance threshold here (between 0 and 1).
    void SetLineNrs(bool val); //!< Set 'true' here to enable line numbers in the output file.
    void SetRegression(bool val); //!< Dummy method for regression (only used for bbrcs).
    //@}
    
    /** @name Others
     *  Other functions.
     */
    //@{
    vector<string>* MineRoot(unsigned int j); //!< Mine fragments rooted at the j-th root node (element type).
    void ReadGsp(FILE* gsp); //!< Read in a gSpan file
    bool AddCompound(string smiles, unsigned int comp_id); //!< Add a compound to the database.
    bool AddActivity(float act, unsigned int comp_id); //!< Add an activity to the database.
    int GetNoRootNodes() {return fm::database->nodelabels.size();} //!< Get number of root nodes (different element types).
    int GetNoCompounds() {return fm::database->trees.size();} //!< Get number of compounds in the database.
    //@}
    
  private:
    void AddChiSqNa(){fm::chisq->na++;fm::chisq->n++;}
    void AddChiSqNi(){fm::chisq->ni++;fm::chisq->n++;}

    bool init_mining_done;
    int comp_runner;
    int comp_no;

    vector<string> r;

};

//#endif
