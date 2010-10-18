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

    unsigned int minfreq;
    int type;
    ChisqConstraint* chisq=NULL;
    // do_backbone missing
    // adjust_ub missing
    bool do_pruning;
    bool aromatic;
    bool refine_singles;
    bool do_output;
    bool bbrc_sep;
    bool regression;


    bool updated;
    // do_yaml missing
    // pvalues missing
    bool gsp_out;
    bool console_out;
    bool db_built;


    bool instance_present;
    int max_hops;


    Database* database=NULL;
    Statistics* statistics=NULL;
    GraphState* graphstate=NULL;


    vector<string>* result=NULL;


    LegOccurrences* legoccurrences=NULL;
    CloseLegOccurrences* closelegoccurrences=NULL; 
    vector<LegOccurrences> candidatelegsoccurrences;
    vector<vector< CloseLegOccurrences> > candidatecloselegsoccs;
    vector<bool> candidatecloselegsoccsused;
    KSConstraint* ks=NULL;


    bool closelegsoccsused;


    // introduced by LAST
    int die;
    bool do_last;
    unsigned int last_hops;



}

#endif
