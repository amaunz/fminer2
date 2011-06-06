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
    unsigned int bbrc_minfreq; // fminer, set
    int bbrc_type;             // fminer, set
    ChisqBbrcConstraint* bbrc_chisq=NULL; // fminer, set (sig).
    bool bbrc_do_backbone; // fminer, set
    bool bbrc_adjust_ub; // fminer, set
    bool bbrc_do_pruning; // fminer, set
    bool bbrc_aromatic; // fminer, set
    bool bbrc_refine_singles; // fminer, set
    bool bbrc_do_output; // fminer, set
    bool bbrc_bbrc_sep; // fminer, set
    bool bbrc_regression; // fminer, set

    // internally controlled by Defaults()
    bool bbrc_updated; // demand
    bool bbrc_do_yaml; // ENV
    bool bbrc_pvalues; // ENV
    bool bbrc_gsp_out; // ENV
    bool bbrc_aromatic_wc; // ENV
    bool bbrc_console_out; // set
    bool bbrc_db_built; // set
    bool bbrc_nr_hits;  // ENV

    // controlled by constructurs & destructor
    bool bbrc_instance_present;

    // controlled by destructor and Reset()
    BbrcDatabase* bbrc_database=NULL;
    BbrcStatistics* bbrc_statistics=NULL;
    BbrcGraphState* bbrc_graphstate=NULL;

    // controlled by Reset()
    vector<string>* bbrc_result=NULL;

    // controlled by destructor & Reset()
    BbrcLegOccurrences* bbrc_legoccurrences=NULL; 
    CloseBbrcLegOccurrences* bbrc_closelegoccurrences=NULL; 
    vector<BbrcLegOccurrences> bbrc_Bbrccandidatelegsoccurrences;
    vector<vector< CloseBbrcLegOccurrences> > bbrc_candidatecloselegsoccs;
    vector<bool> bbrc_candidateBbrccloselegsoccsused;
    KSBbrcConstraint* bbrc_ks=NULL;

    // controlled externally, set on demand
    bool bbrc_Bbrccloselegsoccsused; // demand

}

#endif
