// closeleg.h
// (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010
// Siegfried Nijssen, snijssen@liacs.nl, feb 2004.

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

#ifndef CLOSELEG_H
#define CLOSELEG_H
#include <vector>

#include "misc.h"
#include "legoccurrence.h"

struct BbrcCloseBbrcTuple {
  BbrcEdgeLabel label;
  int from;
  int to;
  friend bool operator< ( BbrcCloseBbrcTuple &a, BbrcCloseBbrcTuple &b ) { return a.from < b.from || ( a.from == b.from && ( a.to < b.to || ( a.to == b.to && a.label < b.label ) ) ); }
  friend bool operator> ( BbrcCloseBbrcTuple &a, BbrcCloseBbrcTuple &b ) { return a.from > b.from || ( a.from == b.from && ( a.to > b.to || ( a.to == b.to && a.label > b.label ) ) ); }
  friend ostream &operator<< ( ostream &stream, BbrcCloseBbrcTuple &tuple ) { 
    stream << (int) tuple.from << " " << tuple.to << " " << (int) tuple.label << endl;
    return stream;
  }
};

struct CloseBbrcLegOccurrence {
  BbrcTid tid;
  BbrcOccurrenceId occurrenceid;

  CloseBbrcLegOccurrence ( BbrcTid tid, BbrcOccurrenceId occurrenceid ): tid ( tid ), occurrenceid ( occurrenceid ) { }
  CloseBbrcLegOccurrence () { }
};

struct CloseBbrcLegOccurrences {
  BbrcFrequency frequency;
  vector<CloseBbrcLegOccurrence> elements;
  CloseBbrcLegOccurrences () : frequency ( 0 ) { }
};

typedef CloseBbrcLegOccurrences *CloseBbrcLegOccurrencesPtr;

struct BbrcCloseBbrcLeg {
  bool copy;
  BbrcCloseBbrcTuple tuple;
  CloseBbrcLegOccurrences occurrences;
  BbrcCloseBbrcLeg (): copy ( true ) { }
};

typedef BbrcCloseBbrcLeg *BbrcCloseBbrcLegPtr;

extern bool Bbrccloselegsoccsused;

class BbrcLeg;
typedef BbrcLeg *BbrcLegPtr;

void BbrcaddCloseExtensions ( vector<BbrcCloseBbrcLegPtr> &targetcloselegs, int number );
void BbrcaddCloseExtensions ( vector<BbrcCloseBbrcLegPtr> &targetcloselegs, vector<BbrcCloseBbrcLegPtr> &sourcecloselegs, BbrcLegOccurrences &sourceoccs );
CloseBbrcLegOccurrencesPtr bbrc_join ( BbrcLegOccurrences &legoccsdata, CloseBbrcLegOccurrences &closelegoccsdata );
CloseBbrcLegOccurrencesPtr bbrc_join ( CloseBbrcLegOccurrences &closelegoccsdata1, CloseBbrcLegOccurrences &closelegoccsdata2 );

#endif
