// globals.h
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

#ifndef GLOBALS_H
#define GLOBALS H

#include "database.h"
#include "constraints.h"

namespace fm {

    unsigned int last_minfreq;
    int last_type;
    ChisqLastConstraint* last_chisq=NULL;
    // do_backbone missing
    // adjust_ub missing
    bool last_do_pruning;
    bool last_aromatic;
    bool last_refine_singles;
    bool last_do_output;
    bool last_bbrc_sep;
    bool last_regression;


    bool last_updated;
    // do_yaml missing
    // pvalues missing
    bool last_gsp_out;
    bool last_console_out;
    bool last_db_built;


    bool last_instance_present;
    int last_max_hops;


    LastDatabase* last_database=NULL;
    LastStatistics* last_statistics=NULL;
    LastGraphState* last_graphstate=NULL;


    vector<string>* last_result=NULL;


    LastLegOccurrences* last_legoccurrences=NULL;
    CloseLastLegOccurrences* last_closelegoccurrences=NULL; 
    vector<LastLegOccurrences> last_Lastcandidatelegsoccurrences;
    vector<vector< CloseLastLegOccurrences> > last_candidatecloselegsoccs;
    vector<bool> last_candidateLastcloselegsoccsused;
    KSLastConstraint* last_ks=NULL;


    bool last_Lastcloselegsoccsused;


    // introduced by LAST
    int last_die;
    bool last_do_last;
    unsigned int last_hops;



}

#endif
