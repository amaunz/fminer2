// closeleg.cpp
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

#include <vector>
#include "misc.h"
#include "closeleg.h"

namespace fm {
    extern unsigned int minfreq;
    extern CloseLegOccurrences* closelegoccurrences;
    extern LegOccurrences* legoccurrences;
    extern vector<vector< CloseLegOccurrences> > candidatecloselegsoccs;
    extern vector<bool> candidatecloselegsoccsused;
    extern bool closelegsoccsused;
}

void addCloseExtensions ( vector<CloseLegPtr> &targetcloselegs, int number ) {
  if ( fm::closelegsoccsused ) {
    for ( int i = 1; i < (int) fm::candidatecloselegsoccs.size (); i++ )
      if ( fm::candidatecloselegsoccsused[i] ) {
        vector<CloseLegOccurrences> &edgelabeloccs = fm::candidatecloselegsoccs[i];
        for ( EdgeLabel j = 0; j < edgelabeloccs.size (); j++ ) {
          if ( edgelabeloccs[j].frequency >= fm::minfreq ) {
            CloseLegPtr closelegptr = new CloseLeg;
            closelegptr->tuple.label = j;
            closelegptr->tuple.to = i;
            closelegptr->tuple.from = number;
            swap ( closelegptr->occurrences, edgelabeloccs[j] );
            targetcloselegs.push_back ( closelegptr );
          }
        }
      }
  }
}

void addCloseExtensions ( vector<CloseLegPtr> &targetcloselegs, vector<CloseLegPtr> &sourcecloselegs, LegOccurrences &sourceoccs ) {
  for ( int i = 0; i < (int) sourcecloselegs.size (); i++ ) {
    CloseLegOccurrencesPtr closelegoccurrencesptr = join ( sourceoccs, sourcecloselegs[i]->occurrences );
    if ( closelegoccurrencesptr ) {
      CloseLegPtr closelegptr = new CloseLeg;
      closelegptr->tuple = sourcecloselegs[i]->tuple;
      swap ( closelegptr->occurrences, *closelegoccurrencesptr );
      targetcloselegs.push_back ( closelegptr );
    }
  }
}

CloseLegOccurrencesPtr join ( LegOccurrences &legoccsdata, CloseLegOccurrences &closelegoccsdata ) {
  Frequency frequency = 0;
  Tid lasttid = NOTID;
  vector<CloseLegOccurrence> &closelegoccs = closelegoccsdata.elements;
  vector<LegOccurrence> &legoccs = legoccsdata.elements;

  fm::closelegoccurrences->elements.resize ( 0 );

  unsigned int legoccssize = legoccs.size (), closelegoccssize = closelegoccs.size ();
  OccurrenceId j = 0, k = 0;
  int comp;

  while ( true ) {
    comp = legoccs[j].occurrenceid - closelegoccs[k].occurrenceid;
    if  ( comp < 0 ) {
      j++;
      if ( j == legoccssize )
        break;
    }
    else {
      if ( comp == 0 ) {
        fm::closelegoccurrences->elements.push_back ( CloseLegOccurrence ( legoccs[j].tid, j ) );
        if ( legoccs[j].tid != lasttid ) {
          lasttid = legoccs[j].tid;
          frequency++;
        }
        j++;
        if ( j == legoccssize )
          break;
      }
      else {
        k++;
        if ( k == closelegoccssize )
          break;
      }
    }
  }

  if ( frequency >= fm::minfreq ) {
    fm::closelegoccurrences->frequency = frequency;
    return fm::closelegoccurrences;
  }
  else
    return NULL;
}

CloseLegOccurrencesPtr join ( CloseLegOccurrences &closelegoccsdata1, CloseLegOccurrences &closelegoccsdata2 ) {
  Frequency frequency = 0;
  Tid lasttid = NOTID;
  vector<CloseLegOccurrence> &closelegoccs1 = closelegoccsdata1.elements,
                             &closelegoccs2 = closelegoccsdata2.elements;

  unsigned int closelegoccs1size = closelegoccs1.size (), closelegoccs2size = closelegoccs2.size ();
  fm::closelegoccurrences->elements.resize ( 0 );
  OccurrenceId j = 0, k = 0;
  int comp;

  while ( true ) {
    comp = closelegoccs1[j].occurrenceid - closelegoccs2[k].occurrenceid;
    if ( comp < 0 ) {
      j++;
      if ( j == closelegoccs1size )
        break;
    }
    else {
      if ( comp == 0 ) {
        fm::closelegoccurrences->elements.push_back ( CloseLegOccurrence ( closelegoccs1[j].tid, closelegoccs1[j].occurrenceid )  );
        if ( closelegoccs1[j].tid != lasttid ) {
          lasttid = closelegoccs1[j].tid;
          frequency++;
        }
        j++;
        if ( j == closelegoccs1size )
          break;
      }
      k++;
      if ( k == closelegoccs2size )
        break;
    }
  }

  if ( frequency >= fm::minfreq ) {
    fm::closelegoccurrences->frequency = frequency;
    return fm::closelegoccurrences;
  }
  else
    return NULL;
}
