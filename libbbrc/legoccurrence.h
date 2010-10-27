// legoccurrence.h
// (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.

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

#ifndef LEGOCCURRENCE_H
#define LEGOCCURRENCE_H
#include <iostream>
#include <vector>

#include "misc.h"

using namespace std;

typedef unsigned int BbrcOccurrenceId;

struct BbrcLegOccurrence {
  BbrcTid tid;
  BbrcOccurrenceId occurrenceid;
  BbrcNodeId tonodeid, fromnodeid;

  BbrcLegOccurrence ( BbrcTid tid, BbrcOccurrenceId occurrenceid, BbrcNodeId tonodeid, BbrcNodeId fromnodeid ): tid ( tid ), occurrenceid ( occurrenceid ), tonodeid ( tonodeid ), fromnodeid ( fromnodeid ) { }
  BbrcLegOccurrence () {}

  friend ostream &operator<< ( ostream &stream, BbrcLegOccurrence &occ );
};

struct BbrcLegOccurrences;
typedef BbrcLegOccurrences *BbrcLegOccurrencesPtr;

struct BbrcLegOccurrences {
  vector<BbrcLegOccurrence> elements;
  BbrcLegOccurrencesPtr parent;
  int number;
  BbrcFrequency selfjoin;
  short unsigned int maxdegree;
  BbrcFrequency frequency;
  BbrcLegOccurrences () : selfjoin ( 0 ), frequency ( 0 ) { }
};

ostream &operator<< ( ostream &stream, vector<BbrcLegOccurrence> &occs );

//extern BbrcLegOccurrences legoccurrences;

// returns the bbrc_join if this bbrc_join is frequent. The returned array may be swapped.
BbrcLegOccurrencesPtr bbrc_join ( BbrcLegOccurrences &legoccsdata1, BbrcNodeId connectingnode, BbrcLegOccurrences &legoccsdata2 );
BbrcLegOccurrencesPtr bbrc_join ( BbrcLegOccurrences &legoccsdata );

extern vector<BbrcLegOccurrences> Bbrccandidatelegsoccurrences; // for each frequent possible edge, the occurrences found, used by bbrc_extend
extern vector<BbrcFrequency> Bbrccandidatelegsfrequencies;

void BbrcinitBbrcLegStatics ();

void bbrc_extend ( BbrcLegOccurrences &legoccurrencesdata ); // fills the global arrays above
void bbrc_extend ( BbrcLegOccurrences &legoccurrencesdata, BbrcEdgeLabel minlabel, BbrcEdgeLabel neglect );

void sanityCheck ( BbrcLegOccurrencesPtr legoccurrencesptr );

#endif
