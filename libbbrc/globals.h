// globals.h
// © 2008 by Andreas Maunz, andreas@maunz.de, jun 2008

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

#ifndef GLOBALS_H
#define GLOBALS H

#include "database.h"
#include "constraints.h"

namespace fm {

    unsigned int minfreq;
    int type;
    bool do_backbone;
    bool updated;
    bool adjust_ub;
    bool do_pruning;
    bool instance_present;
    bool console_out;
    bool aromatic;
    bool refine_singles;
    bool do_output;
    bool do_yaml;
    bool gsp_out;
    bool pvalues;
    bool bbrc_sep;
    bool most_specific_trees_only;
    bool line_nrs;
    bool regression;

    Database* database=NULL;
    Statistics* statistics=NULL;
    ChisqConstraint* chisq=NULL;
    KSConstraint* ks=NULL;
    GraphState* graphstate=NULL;
    CloseLegOccurrences* closelegoccurrences=NULL; 
    LegOccurrences* legoccurrences=NULL;

    vector<string>* result=NULL;
    vector<LegOccurrences> candidatelegsoccurrences;
    vector<vector< CloseLegOccurrences> > candidatecloselegsoccs;
    vector<bool> candidatecloselegsoccsused;

    bool closelegsoccsused;

}

#endif
