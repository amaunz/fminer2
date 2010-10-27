// closeleg.h
// (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010
// Siegfried Nijssen, snijssen@liacs.nl, feb 2004.

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

#ifndef CLOSELEG_H
#define CLOSELEG_H
#include <vector>

#include "misc.h"
#include "legoccurrence.h"

struct LastCloseLastTuple {
  LastEdgeLabel label;
  int from;
  int to;
  friend bool operator< ( LastCloseLastTuple &a, LastCloseLastTuple &b ) { return a.from < b.from || ( a.from == b.from && ( a.to < b.to || ( a.to == b.to && a.label < b.label ) ) ); }
  friend bool operator> ( LastCloseLastTuple &a, LastCloseLastTuple &b ) { return a.from > b.from || ( a.from == b.from && ( a.to > b.to || ( a.to == b.to && a.label > b.label ) ) ); }
  friend ostream &operator<< ( ostream &stream, LastCloseLastTuple &tuple ) { 
    stream << (int) tuple.from << " " << tuple.to << " " << (int) tuple.label << endl;
    return stream;
  }
};

struct CloseLastLastLegOccurrence {
  LastTid tid;
  LastOccurrenceId occurrenceid;

  CloseLastLastLegOccurrence ( LastTid tid, LastOccurrenceId occurrenceid ): tid ( tid ), occurrenceid ( occurrenceid ) { }
  CloseLastLastLegOccurrence () { }
};

struct CloseLastLastLegOccurrences {
  LastFrequency frequency;
  vector<CloseLastLastLegOccurrence> elements;
  CloseLastLastLegOccurrences () : frequency ( 0 ) { }
};

typedef CloseLastLastLegOccurrences *CloseLastLastLegOccurrencesPtr;

struct LastCloseLastLeg {
  bool copy;
  LastCloseLastTuple tuple;
  CloseLastLastLegOccurrences occurrences;
  LastCloseLastLeg (): copy ( true ) { }
};

typedef LastCloseLastLeg *LastCloseLastLegPtr;

extern bool Lastcloselegsoccsused;

class LastLeg;
typedef LastLeg *LastLegPtr;

void LastaddCloseExtensions ( vector<LastCloseLastLegPtr> &targetcloselegs, int number );
void LastaddCloseExtensions ( vector<LastCloseLastLegPtr> &targetcloselegs, vector<LastCloseLastLegPtr> &sourcecloselegs, LastLastLegOccurrences &sourceoccs );
CloseLastLastLegOccurrencesPtr join ( LastLastLegOccurrences &legoccsdata, CloseLastLastLegOccurrences &closelegoccsdata );
CloseLastLastLegOccurrencesPtr join ( CloseLastLastLegOccurrences &closelegoccsdata1, CloseLastLastLegOccurrences &closelegoccsdata2 );

#endif
