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

struct BbrcBbrcLegOccurrence {
  BbrcTid tid;
  BbrcOccurrenceId occurrenceid;
  BbrcNodeId tonodeid, fromnodeid;

  BbrcBbrcLegOccurrence ( BbrcTid tid, BbrcOccurrenceId occurrenceid, BbrcNodeId tonodeid, BbrcNodeId fromnodeid ): tid ( tid ), occurrenceid ( occurrenceid ), tonodeid ( tonodeid ), fromnodeid ( fromnodeid ) { }
  BbrcBbrcLegOccurrence () {}

  friend ostream &operator<< ( ostream &stream, BbrcBbrcLegOccurrence &occ );
};

struct BbrcBbrcLegOccurrences;
typedef BbrcBbrcLegOccurrences *BbrcBbrcLegOccurrencesPtr;

struct BbrcBbrcLegOccurrences {
  vector<BbrcBbrcLegOccurrence> elements;
  BbrcBbrcLegOccurrencesPtr parent;
  int number;
  BbrcFrequency selfjoin;
  short unsigned int maxdegree;
  BbrcFrequency frequency;
  BbrcBbrcLegOccurrences () : selfjoin ( 0 ), frequency ( 0 ) { }
};

ostream &operator<< ( ostream &stream, vector<BbrcBbrcLegOccurrence> &occs );

//extern BbrcBbrcLegOccurrences legoccurrences;

// returns the bbrc_join if this bbrc_join is frequent. The returned array may be swapped.
BbrcBbrcLegOccurrencesPtr bbrc_join ( BbrcBbrcLegOccurrences &legoccsdata1, BbrcNodeId connectingnode, BbrcBbrcLegOccurrences &legoccsdata2 );
BbrcBbrcLegOccurrencesPtr bbrc_join ( BbrcBbrcLegOccurrences &legoccsdata );

extern vector<BbrcBbrcLegOccurrences> Bbrccandidatelegsoccurrences; // for each frequent possible edge, the occurrences found, used by bbrc_extend
extern vector<BbrcFrequency> Bbrccandidatelegsfrequencies;

void BbrcinitBbrcLegStatics ();

void bbrc_extend ( BbrcBbrcLegOccurrences &legoccurrencesdata ); // fills the global arrays above
void bbrc_extend ( BbrcBbrcLegOccurrences &legoccurrencesdata, BbrcEdgeLabel minlabel, BbrcEdgeLabel neglect );

void sanityCheck ( BbrcBbrcLegOccurrencesPtr legoccurrencesptr );

#endif
