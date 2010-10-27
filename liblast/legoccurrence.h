// legoccurrence.h
// (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.

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

#ifndef LEGOCCURRENCE_H
#define LEGOCCURRENCE_H
#include <iostream>
#include <vector>

#include "misc.h"

using namespace std;

typedef unsigned int LastOccurrenceId;

struct LastLegOccurrence {
  LastTid tid;
  LastOccurrenceId occurrenceid;
  LastNodeId tonodeid, fromnodeid;

  LastLegOccurrence ( LastTid tid, LastOccurrenceId occurrenceid, LastNodeId tonodeid, LastNodeId fromnodeid ): tid ( tid ), occurrenceid ( occurrenceid ), tonodeid ( tonodeid ), fromnodeid ( fromnodeid ) { }
  LastLegOccurrence () {}

  friend ostream &operator<< ( ostream &stream, LastLegOccurrence &occ );
};

struct LastLegOccurrences;
typedef LastLegOccurrences *LastLegOccurrencesPtr;

struct LastLegOccurrences {
  vector<LastLegOccurrence> elements;
  LastLegOccurrencesPtr parent;
  int number;
  LastFrequency selfjoin;
  short unsigned int maxdegree;
  LastFrequency frequency;
  LastLegOccurrences () : selfjoin ( 0 ), frequency ( 0 ) { }
};

ostream &operator<< ( ostream &stream, vector<LastLegOccurrence> &occs );

//extern LastLegOccurrences legoccurrences;

// returns the bbrc_join if this bbrc_join is frequent. The returned array may be swapped.
LastLegOccurrencesPtr bbrc_join ( LastLegOccurrences &legoccsdata1, LastNodeId connectingnode, LastLegOccurrences &legoccsdata2 );
LastLegOccurrencesPtr bbrc_join ( LastLegOccurrences &legoccsdata );

extern vector<LastLegOccurrences> Lastcandidatelegsoccurrences; // for each frequent possible edge, the occurrences found, used by bbrc_extend
extern vector<LastFrequency> Lastcandidatelegsfrequencies;

void LastinitLastLegStatics ();

void bbrc_extend ( LastLegOccurrences &legoccurrencesdata ); // fills the global arrays above
void bbrc_extend ( LastLegOccurrences &legoccurrencesdata, LastEdgeLabel minlabel, LastEdgeLabel neglect );

void sanityCheck ( LastLegOccurrencesPtr legoccurrencesptr );

#endif
