// closeleg.cpp
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

#include <vector>
#include "misc.h"
#include "closeleg.h"

namespace fm {
    extern float last_minfreq;
    extern CloseLastLegOccurrences* last_closelegoccurrences;
    extern LastLegOccurrences* last_legoccurrences;
    extern vector<vector< CloseLastLegOccurrences> > last_candidatecloselegsoccs;
    extern vector<bool> last_candidateLastcloselegsoccsused;
    extern bool last_Lastcloselegsoccsused;
}

void LastaddCloseExtensions ( vector<LastCloseLastLegPtr> &targetcloselegs, int number ) {
  if ( fm::last_Lastcloselegsoccsused ) {
    for ( int i = 1; i < (int) fm::last_candidatecloselegsoccs.size (); i++ )
      if ( fm::last_candidateLastcloselegsoccsused[i] ) {
        vector<CloseLastLegOccurrences> &edgelabeloccs = fm::last_candidatecloselegsoccs[i];
        for ( LastEdgeLabel j = 0; j < edgelabeloccs.size (); j++ ) {
          if ( edgelabeloccs[j].frequency >= fm::last_minfreq ) {
            LastCloseLastLegPtr closelegptr = new LastCloseLastLeg;
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

void LastaddCloseExtensions ( vector<LastCloseLastLegPtr> &targetcloselegs, vector<LastCloseLastLegPtr> &sourcecloselegs, LastLegOccurrences &sourceoccs ) {
  for ( int i = 0; i < (int) sourcecloselegs.size (); i++ ) {
    CloseLastLegOccurrencesPtr closelegoccurrencesptr = bbrc_join ( sourceoccs, sourcecloselegs[i]->occurrences );
    if ( closelegoccurrencesptr ) {
      LastCloseLastLegPtr closelegptr = new LastCloseLastLeg;
      closelegptr->tuple = sourcecloselegs[i]->tuple;
      swap ( closelegptr->occurrences, *closelegoccurrencesptr );
      targetcloselegs.push_back ( closelegptr );
    }
  }
}

CloseLastLegOccurrencesPtr bbrc_join ( LastLegOccurrences &legoccsdata, CloseLastLegOccurrences &closelegoccsdata ) {
  LastFrequency frequency = 0;
  LastTid lasttid = NOTID;
  vector<CloseLastLegOccurrence> &closelegoccs = closelegoccsdata.elements;
  vector<LastLegOccurrence> &legoccs = legoccsdata.elements;

  fm::last_closelegoccurrences->elements.resize ( 0 );

  unsigned int legoccssize = legoccs.size (), closelegoccssize = closelegoccs.size ();
  LastOccurrenceId j = 0, k = 0;
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
        fm::last_closelegoccurrences->elements.push_back ( CloseLastLegOccurrence ( legoccs[j].tid, j ) );
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

  if ( frequency >= fm::last_minfreq ) {
    fm::last_closelegoccurrences->frequency = frequency;
    return fm::last_closelegoccurrences;
  }
  else
    return NULL;
}

CloseLastLegOccurrencesPtr bbrc_join ( CloseLastLegOccurrences &closelegoccsdata1, CloseLastLegOccurrences &closelegoccsdata2 ) {
  LastFrequency frequency = 0;
  LastTid lasttid = NOTID;
  vector<CloseLastLegOccurrence> &closelegoccs1 = closelegoccsdata1.elements,
                             &closelegoccs2 = closelegoccsdata2.elements;

  unsigned int closelegoccs1size = closelegoccs1.size (), closelegoccs2size = closelegoccs2.size ();
  fm::last_closelegoccurrences->elements.resize ( 0 );
  LastOccurrenceId j = 0, k = 0;
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
        fm::last_closelegoccurrences->elements.push_back ( CloseLastLegOccurrence ( closelegoccs1[j].tid, closelegoccs1[j].occurrenceid )  );
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

  if ( frequency >= fm::last_minfreq ) {
    fm::last_closelegoccurrences->frequency = frequency;
    return fm::last_closelegoccurrences;
  }
  else
    return NULL;
}
