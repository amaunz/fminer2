// closeleg.cpp
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

#include <vector>
#include "misc.h"
#include "closeleg.h"
#include "database.h"

namespace fm {
    extern BbrcDatabase* bbrc_database;
    extern float bbrc_minfreq;
    extern CloseBbrcLegOccurrences* bbrc_closelegoccurrences;
    extern BbrcLegOccurrences* bbrc_legoccurrences;
    extern vector<vector< CloseBbrcLegOccurrences> > bbrc_candidatecloselegsoccs;
    extern vector<bool> bbrc_candidateBbrccloselegsoccsused;
    extern bool bbrc_Bbrccloselegsoccsused;
}

void BbrcaddCloseExtensions ( vector<BbrcCloseBbrcLegPtr> &targetcloselegs, int number ) {
  if ( fm::bbrc_Bbrccloselegsoccsused ) {
    for ( int i = 1; i < (int) fm::bbrc_candidatecloselegsoccs.size (); i++ )
      if ( fm::bbrc_candidateBbrccloselegsoccsused[i] ) {
        vector<CloseBbrcLegOccurrences> &edgelabeloccs = fm::bbrc_candidatecloselegsoccs[i];
        for ( BbrcEdgeLabel j = 0; j < edgelabeloccs.size (); j++ ) {
          if ( edgelabeloccs[j].frequency >= fm::bbrc_minfreq ) {
            BbrcCloseBbrcLegPtr closelegptr = new BbrcCloseBbrcLeg;
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

void BbrcaddCloseExtensions ( vector<BbrcCloseBbrcLegPtr> &targetcloselegs, vector<BbrcCloseBbrcLegPtr> &sourcecloselegs, BbrcLegOccurrences &sourceoccs ) {
  for ( int i = 0; i < (int) sourcecloselegs.size (); i++ ) {
    CloseBbrcLegOccurrencesPtr closelegoccurrencesptr = bbrc_join ( sourceoccs, sourcecloselegs[i]->occurrences );
    if ( closelegoccurrencesptr ) {
      BbrcCloseBbrcLegPtr closelegptr = new BbrcCloseBbrcLeg;
      closelegptr->tuple = sourcecloselegs[i]->tuple;
      swap ( closelegptr->occurrences, *closelegoccurrencesptr );
      targetcloselegs.push_back ( closelegptr );
    }
  }
}

CloseBbrcLegOccurrencesPtr bbrc_join ( BbrcLegOccurrences &legoccsdata, CloseBbrcLegOccurrences &closelegoccsdata ) {
  BbrcFrequency frequency = 0;
  BbrcTid lasttid = NOTID;
  vector<CloseBbrcLegOccurrence> &closelegoccs = closelegoccsdata.elements;
  vector<BbrcLegOccurrence> &legoccs = legoccsdata.elements;

  fm::bbrc_closelegoccurrences->elements.resize ( 0 );

  unsigned int legoccssize = legoccs.size (), closelegoccssize = closelegoccs.size ();
  BbrcOccurrenceId j = 0, k = 0;
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
        fm::bbrc_closelegoccurrences->elements.push_back ( CloseBbrcLegOccurrence ( legoccs[j].tid, j ) );
        if ( legoccs[j].tid != lasttid ) {
          lasttid = legoccs[j].tid;
          //frequency++;
          frequency += 
          fm::bbrc_database->trees[legoccs[j].tid]->weight;
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

  if ( frequency >= fm::bbrc_minfreq ) {
    fm::bbrc_closelegoccurrences->frequency = frequency;
    return fm::bbrc_closelegoccurrences;
  }
  else
    return NULL;
}

CloseBbrcLegOccurrencesPtr bbrc_join ( CloseBbrcLegOccurrences &closelegoccsdata1, CloseBbrcLegOccurrences &closelegoccsdata2 ) {
  BbrcFrequency frequency = 0;
  BbrcTid lasttid = NOTID;
  vector<CloseBbrcLegOccurrence> &closelegoccs1 = closelegoccsdata1.elements,
                             &closelegoccs2 = closelegoccsdata2.elements;

  unsigned int closelegoccs1size = closelegoccs1.size (), closelegoccs2size = closelegoccs2.size ();
  fm::bbrc_closelegoccurrences->elements.resize ( 0 );
  BbrcOccurrenceId j = 0, k = 0;
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
        fm::bbrc_closelegoccurrences->elements.push_back ( CloseBbrcLegOccurrence ( closelegoccs1[j].tid, closelegoccs1[j].occurrenceid )  );
        if ( closelegoccs1[j].tid != lasttid ) {
          lasttid = closelegoccs1[j].tid;
          // frequency++;
          frequency += 
          fm::bbrc_database->trees[closelegoccs1[j].tid]->weight;
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

  if ( frequency >= fm::bbrc_minfreq ) {
    fm::bbrc_closelegoccurrences->frequency = frequency;
    return fm::bbrc_closelegoccurrences;
  }
  else
    return NULL;
}
