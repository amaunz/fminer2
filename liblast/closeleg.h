// closeleg.h
// © 2008 by Andreas Maunz, andreas@maunz.de, jul 2008
// Siegfried Nijssen, snijssen@liacs.nl, feb 2004.

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

#ifndef CLOSELEG_H
#define CLOSELEG_H
#include <vector>

#include "misc.h"
#include "legoccurrence.h"

struct CloseTuple {
  EdgeLabel label;
  int from;
  int to;
  friend bool operator< ( CloseTuple &a, CloseTuple &b ) { return a.from < b.from || ( a.from == b.from && ( a.to < b.to || ( a.to == b.to && a.label < b.label ) ) ); }
  friend bool operator> ( CloseTuple &a, CloseTuple &b ) { return a.from > b.from || ( a.from == b.from && ( a.to > b.to || ( a.to == b.to && a.label > b.label ) ) ); }
  friend ostream &operator<< ( ostream &stream, CloseTuple &tuple ) { 
    stream << (int) tuple.from << " " << tuple.to << " " << (int) tuple.label << endl;
    return stream;
  }
};

struct CloseLegOccurrence {
  Tid tid;
  OccurrenceId occurrenceid;

  CloseLegOccurrence ( Tid tid, OccurrenceId occurrenceid ): tid ( tid ), occurrenceid ( occurrenceid ) { }
  CloseLegOccurrence () { }
};

struct CloseLegOccurrences {
  Frequency frequency;
  vector<CloseLegOccurrence> elements;
  CloseLegOccurrences () : frequency ( 0 ) { }
};

typedef CloseLegOccurrences *CloseLegOccurrencesPtr;

struct CloseLeg {
  bool copy;
  CloseTuple tuple;
  CloseLegOccurrences occurrences;
  CloseLeg (): copy ( true ) { }
};

typedef CloseLeg *CloseLegPtr;

extern bool closelegsoccsused;

class Leg;
typedef Leg *LegPtr;

void addCloseExtensions ( vector<CloseLegPtr> &targetcloselegs, int number );
void addCloseExtensions ( vector<CloseLegPtr> &targetcloselegs, vector<CloseLegPtr> &sourcecloselegs, LegOccurrences &sourceoccs );
CloseLegOccurrencesPtr join ( LegOccurrences &legoccsdata, CloseLegOccurrences &closelegoccsdata );
CloseLegOccurrencesPtr join ( CloseLegOccurrences &closelegoccsdata1, CloseLegOccurrences &closelegoccsdata2 );

#endif
