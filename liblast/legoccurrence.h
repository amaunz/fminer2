// legoccurrence.h
// © 2008 by Andreas Maunz, andreas@maunz.de, jul 2008
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.

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

#ifndef LEGOCCURRENCE_H
#define LEGOCCURRENCE_H
#include <iostream>
#include <vector>

#include "misc.h"

using namespace std;

typedef unsigned int OccurrenceId;

struct LegOccurrence {
  Tid tid;
  OccurrenceId occurrenceid;
  NodeId tonodeid, fromnodeid;

  LegOccurrence ( Tid tid, OccurrenceId occurrenceid, NodeId tonodeid, NodeId fromnodeid ): tid ( tid ), occurrenceid ( occurrenceid ), tonodeid ( tonodeid ), fromnodeid ( fromnodeid ) { }
  LegOccurrence () {}

  friend ostream &operator<< ( ostream &stream, LegOccurrence &occ );
};

struct LegOccurrences;
typedef LegOccurrences *LegOccurrencesPtr;

struct LegOccurrences {
  vector<LegOccurrence> elements;
  LegOccurrencesPtr parent;
  int number;
  Frequency selfjoin;
  short unsigned int maxdegree;
  Frequency frequency;
  LegOccurrences () : selfjoin ( 0 ), frequency ( 0 ) { }
};

ostream &operator<< ( ostream &stream, vector<LegOccurrence> &occs );

//extern LegOccurrences legoccurrences;

// returns the join if this join is frequent. The returned array may be swapped.
LegOccurrencesPtr join ( LegOccurrences &legoccsdata1, NodeId connectingnode, LegOccurrences &legoccsdata2 );
LegOccurrencesPtr join ( LegOccurrences &legoccsdata );

extern vector<LegOccurrences> candidatelegsoccurrences; // for each frequent possible edge, the occurrences found, used by extend
extern vector<Frequency> candidatelegsfrequencies;

void initLegStatics ();

void extend ( LegOccurrences &legoccurrencesdata ); // fills the global arrays above
void extend ( LegOccurrences &legoccurrencesdata, EdgeLabel minlabel, EdgeLabel neglect );

void sanityCheck ( LegOccurrencesPtr legoccurrencesptr );

#endif
