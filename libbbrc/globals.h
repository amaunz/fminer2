// globals.h
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

#ifndef GLOBALS_H
#define GLOBALS H

#include "database.h"
#include "constraints.h"

namespace fm {
    // switched by fminer binary
    unsigned int minfreq; // fminer, set
    int type;             // fminer, set
    ChisqConstraint* chisq=NULL; // fminer, set (sig).
    bool do_backbone; // fminer, set
    bool adjust_ub; // fminer, set
    bool do_pruning; // fminer, set
    bool aromatic; // fminer, set
    bool refine_singles; // fminer, set
    bool do_output; // fminer, set
    bool bbrc_sep; // fminer, set
    bool regression; // fminer, set

    // internally controlled by Defaults()
    bool updated; // demand
    bool do_yaml; // ENV
    bool pvalues; // ENV
    bool gsp_out; // ENV
    bool no_aromatic; // ENV
    bool console_out; // set

    // controlled by constructurs & destructor
    bool instance_present;

    // controlled by destructor and Reset()
    Database* database=NULL;
    Statistics* statistics=NULL;
    GraphState* graphstate=NULL;

    // controlled by Reset()
    vector<string>* result=NULL;

    // controlled by destructor & Reset()
    LegOccurrences* legoccurrences=NULL; 
    CloseLegOccurrences* closelegoccurrences=NULL; 
    vector<LegOccurrences> candidatelegsoccurrences;
    vector<vector< CloseLegOccurrences> > candidatecloselegsoccs;
    vector<bool> candidatecloselegsoccsused;
    KSConstraint* ks=NULL;

    // controlled externally, set on demand
    bool closelegsoccsused; // demand

}

#endif
