// legoccurrence.cpp
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

#include "legoccurrence.h"
#include "closeleg.h"
#include "database.h"
#include "graphstate.h"

namespace fm {
    extern BbrcDatabase* database;
    extern BbrcGraphState* graphstate;
    extern BbrcLegOccurrences* legoccurrences;
    extern unsigned int minfreq;
    extern vector<BbrcLegOccurrences> Bbrccandidatelegsoccurrences; 
    extern vector<vector< CloseBbrcLegOccurrences> > candidatecloselegsoccs;
    extern vector<bool> candidateBbrccloselegsoccsused;
    extern bool Bbrccloselegsoccsused;
}

void BbrcinitBbrcLegStatics () {
  fm::candidatecloselegsoccs.reserve ( 200 ); // should be larger than the largest structure that contains a cycle
  fm::Bbrccandidatelegsoccurrences.resize ( fm::database->frequentBbrcEdgeLabelSize () );
}


/*
ostream &operator<< ( ostream &stream, BbrcLegOccurrence &occ ) {
  stream << "[" << occ.tid << "(" << fm::database->trees[occ.tid]->activity  << ")" << "," << occ.occurrenceid << "," << occ.tonodeid << "," << occ.fromnodeid << "]";
  return stream;
}
*/

ostream &operator<< ( ostream &stream, vector<BbrcLegOccurrence> &occs ) {
  BbrcTid lasttid = NOTID;
  BbrcFrequency frequency = 0;
  for ( int i = 0; i < (int) occs.size (); i++ ) {
    //stream << occs[i];
    if ( occs[i].tid != lasttid ) {
      stream << occs[i].tid << " ";
      lasttid = occs[i].tid;
      frequency++;
    }
  }
  stream << endl << " (" << frequency << ")" << endl;

  return stream;
}

// This function is on the critical path. Its efficiency is MOST important.
BbrcLegOccurrencesPtr bbrc_join ( BbrcLegOccurrences &legoccsdata1, BbrcNodeId connectingnode, BbrcLegOccurrences &legoccsdata2 ) {
  if ( fm::graphstate->getNodeDegree ( connectingnode ) == fm::graphstate->getNodeMaxDegree ( connectingnode ) ) 
    return NULL;

  BbrcFrequency frequency = 0;
  BbrcTid lasttid = NOTID;
  vector<BbrcLegOccurrence> &legoccs1 = legoccsdata1.elements, &legoccs2 = legoccsdata2.elements;
  fm::legoccurrences->elements.resize ( 0 );
  fm::legoccurrences->maxdegree = 0;
  fm::legoccurrences->selfjoin = 0;
  //fm::legoccurrences->elements.reserve ( legoccs1.size () * 2 ); // increased memory usage, and speed!
  BbrcOccurrenceId j = 0, k = 0, l, m;
  unsigned int legoccs1size = legoccs1.size (), legoccs2size = legoccs2.size (); // this increases speed CONSIDERABLY!
  BbrcTid lastself = NOTID;

  do {
    while ( j < legoccs1size && legoccs1[j].occurrenceid < legoccs2[k].occurrenceid ) {
      j++;
    }
    if ( j < legoccs1size ) {
      BbrcLegOccurrence &jlegocc = legoccs1[j];
      while ( k < legoccs2size && legoccs2[k].occurrenceid < jlegocc.occurrenceid ) {
        k++;
      }
      if ( k < legoccs2size ) {
        if ( legoccs2[k].occurrenceid == jlegocc.occurrenceid ) {
          m = j;
          do {
            j++;
          }
          while ( j < legoccs1size && legoccs1[j].occurrenceid == jlegocc.occurrenceid );
          l = k;
          do {
            k++;
          }
          while ( k < legoccs2size && legoccs2[k].occurrenceid == jlegocc.occurrenceid );
    	  bool add = false;
          for ( BbrcOccurrenceId m2 = m; m2 < j; m2++ ) {
            int d = 0;
            for ( BbrcOccurrenceId l2 = l; l2 < k; l2++ ) {
	      BbrcNodeId tonodeid = legoccs2[l2].tonodeid;
              if ( legoccs1[m2].tonodeid !=  tonodeid ) {
                fm::legoccurrences->elements.push_back ( BbrcLegOccurrence ( jlegocc.tid, m2, tonodeid, legoccs2[l2].fromnodeid ) );
                Bbrcsetmax ( fm::legoccurrences->maxdegree, fm::database->trees[jlegocc.tid]->nodes[tonodeid].edges.size () );
        		add = true;
        		d++;
              }
            }
	    if ( d > 1 && jlegocc.tid != lastself ) {
	      lastself = jlegocc.tid;
	      fm::legoccurrences->selfjoin++;
	    }
	  }
	  	  
	  if ( jlegocc.tid != lasttid && add ) {
        lasttid = jlegocc.tid;
	    frequency++;
	  }

          if ( k == legoccs2size )
            break;
        }
      }
      else
        break;
    }
    else
      break;
  }
  while ( true );

  if ( frequency >= fm::minfreq ) {
    fm::legoccurrences->parent = &legoccsdata1;
    fm::legoccurrences->number = legoccsdata1.number + 1;
    fm::legoccurrences->frequency = frequency;
    return fm::legoccurrences;
  }
  else
    return NULL;
}

BbrcLegOccurrencesPtr bbrc_join ( BbrcLegOccurrences &legoccsdata ) {
  if ( legoccsdata.selfjoin < fm::minfreq ) 
    return NULL;
  fm::legoccurrences->elements.resize ( 0 );
  vector<BbrcLegOccurrence> &legoccs = legoccsdata.elements;
  fm::legoccurrences->maxdegree = 0;
  fm::legoccurrences->selfjoin = 0;
  BbrcTid lastself = NOTID;

  BbrcOccurrenceId j = 0, k, l, m;
  do {
    k = j;
    BbrcLegOccurrence &legocc = legoccs[k];
    do {
      j++;
    }
    while ( j < legoccs.size () &&
            legoccs[j].occurrenceid == legocc.occurrenceid );
    for ( l = k; l < j; l++ )
      for ( m = k; m < j; m++ )
        if ( l != m ) {
          fm::legoccurrences->elements.push_back ( BbrcLegOccurrence ( legocc.tid, l, legoccs[m].tonodeid, legoccs[m].fromnodeid ) );
          Bbrcsetmax ( fm::legoccurrences->maxdegree, fm::database->trees[legocc.tid]->nodes[legoccs[m].tonodeid].edges.size () );
        }
    if ( ( j - k > 2 ) && legocc.tid != lastself ) {
      lastself = legocc.tid;
      fm::legoccurrences->selfjoin++;
    }
  }
  while ( j < legoccs.size () );

    // no need to check that we are frequent, we must be frequent
  fm::legoccurrences->parent = &legoccsdata;
  fm::legoccurrences->number = legoccsdata.number + 1;
  fm::legoccurrences->frequency = legoccsdata.selfjoin; 
    // we compute the self-join frequency exactly while building the
    // previous list. It is therefore not necessary to recompute it.
  return fm::legoccurrences;
}

inline int nocycle ( BbrcDatabaseTreePtr tree, BbrcDatabaseTreeNode &node, BbrcNodeId tonode, BbrcOccurrenceId occurrenceid, BbrcLegOccurrencesPtr legoccurrencesdataptr ) {
  if ( !tree->nodes[tonode].incycle )
    return 0;
  if ( !node.incycle )
    return 0;
  while ( legoccurrencesdataptr ) {
    if ( legoccurrencesdataptr->elements[occurrenceid].tonodeid == tonode ) {
      return legoccurrencesdataptr->number;
    }
    occurrenceid = legoccurrencesdataptr->elements[occurrenceid].occurrenceid;
    legoccurrencesdataptr = legoccurrencesdataptr->parent;
  }
  return 0;
}

void candidateBbrcCloseBbrcLegsAllocate ( int number, int maxnumber ) {
  if ( !fm::Bbrccloselegsoccsused ) {
    int oldsize = fm::candidatecloselegsoccs.size ();
    fm::candidatecloselegsoccs.resize ( maxnumber );
    for ( int k = oldsize; k < (int) fm::candidatecloselegsoccs.size (); k++ ) {
      fm::candidatecloselegsoccs[k].resize ( fm::database->frequentBbrcEdgeLabelSize () );
    }
    fm::candidateBbrccloselegsoccsused.resize ( 0 );
    fm::candidateBbrccloselegsoccsused.resize ( maxnumber, false );
    fm::Bbrccloselegsoccsused = true;
  }
  if ( !fm::candidateBbrccloselegsoccsused[number] ) {
    fm::candidateBbrccloselegsoccsused[number] = true;
    vector<CloseBbrcLegOccurrences> &candidateedgelabeloccs = fm::candidatecloselegsoccs[number];
    for ( int k = 0; k < (int) candidateedgelabeloccs.size (); k++ ) {
      candidateedgelabeloccs[k].elements.resize ( 0 );
      candidateedgelabeloccs[k].frequency = 0;
    }
  }
}






void bbrc_extend ( BbrcLegOccurrences &legoccurrencesdata ) {
  // we're trying hard to avoid repeated destructor/constructor calls for complex types like vectors.
  // better reuse previously allocated memory, if possible!
  
  

  vector<BbrcLegOccurrence> &legoccurrences = legoccurrencesdata.elements;   ///////////////////////////////////////AM : BUG!!!



  BbrcTid lastself[fm::Bbrccandidatelegsoccurrences.size ()];

  for ( int i = 0; i < (int) fm::Bbrccandidatelegsoccurrences.size (); i++ ) {
    fm::Bbrccandidatelegsoccurrences[i].elements.resize ( 0 );
    //fm::Bbrccandidatelegsoccurrences[i].elements.reserve ( legoccurrences.size () ); // increases memory usage, but also speed!
    fm::Bbrccandidatelegsoccurrences[i].parent = &legoccurrencesdata;
    fm::Bbrccandidatelegsoccurrences[i].number = legoccurrencesdata.number + 1;
    fm::Bbrccandidatelegsoccurrences[i].maxdegree = 0;
    fm::Bbrccandidatelegsoccurrences[i].frequency = 0;
    fm::Bbrccandidatelegsoccurrences[i].selfjoin = 0;
    lastself[i] = NOTID;
  }

  fm::Bbrccloselegsoccsused = false; // we are lazy with the initialization of close leg arrays, as we may not need them at all in
                             // many cases

  for ( BbrcOccurrenceId i = 0; i < legoccurrences.size (); i++ ) {
    BbrcLegOccurrence &legocc = legoccurrences[i];
    BbrcDatabaseTreePtr tree = fm::database->trees[legocc.tid];
    BbrcDatabaseTreeNode &node = tree->nodes[legocc.tonodeid];
    for ( int j = 0; j < node.edges.size (); j++ ) {
      if ( node.edges[j].tonode != legocc.fromnodeid ) {
      	BbrcEdgeLabel edgelabel = node.edges[j].edgelabel;

        int number = nocycle ( tree, node, node.edges[j].tonode, i, &legoccurrencesdata );

        if ( number == 0 ) {
          vector<BbrcLegOccurrence> &candidatelegsoccs = fm::Bbrccandidatelegsoccurrences[edgelabel].elements;
          if ( candidatelegsoccs.empty () )  fm::Bbrccandidatelegsoccurrences[edgelabel].frequency++;
          else {

	            if ( candidatelegsoccs.back ().tid != legocc.tid )
        	        fm::Bbrccandidatelegsoccurrences[edgelabel].frequency++;

	            if ( candidatelegsoccs.back ().occurrenceid == i &&
	                lastself[edgelabel] != legocc.tid ) {
                    lastself[edgelabel] = legocc.tid;
	                fm::Bbrccandidatelegsoccurrences[edgelabel].selfjoin++;
	            }

          }
          candidatelegsoccs.push_back ( BbrcLegOccurrence ( legocc.tid, i, node.edges[j].tonode, legocc.tonodeid ) );
          Bbrcsetmax ( fm::Bbrccandidatelegsoccurrences[edgelabel].maxdegree, fm::database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
        }

        else if ( number - 1 != fm::graphstate->nodes.back().edges[0].tonode ) {
            candidateBbrcCloseBbrcLegsAllocate ( number, legoccurrencesdata.number + 1 );
            vector<CloseBbrcLegOccurrence> &candidatelegsoccs = fm::candidatecloselegsoccs[number][edgelabel].elements;
            if ( !candidatelegsoccs.size () || candidatelegsoccs.back ().tid != legocc.tid )
	            fm::candidatecloselegsoccs[number][edgelabel].frequency++;
            candidatelegsoccs.push_back ( CloseBbrcLegOccurrence ( legocc.tid, i ) );
            Bbrcsetmax ( fm::Bbrccandidatelegsoccurrences[edgelabel].maxdegree, fm::database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
        }

      }
    }
  }
}







void bbrc_extend ( BbrcLegOccurrences &legoccurrencesdata, BbrcEdgeLabel minlabel, BbrcEdgeLabel neglect ) {


  // we're trying hard to avoid repeated destructor/constructor calls for complex types like vectors.
  // better reuse previously allocated memory, if possible!




  vector<BbrcLegOccurrence> &legoccurrences = legoccurrencesdata.elements;  ///////////////////////////////////////AM : BUG!!!




  int lastself[fm::Bbrccandidatelegsoccurrences.size ()];
  
  for ( int i = 0; i < (int) fm::Bbrccandidatelegsoccurrences.size (); i++ ) {
    fm::Bbrccandidatelegsoccurrences[i].elements.resize ( 0 );
    fm::Bbrccandidatelegsoccurrences[i].parent = &legoccurrencesdata;
    fm::Bbrccandidatelegsoccurrences[i].number = legoccurrencesdata.number + 1;
    fm::Bbrccandidatelegsoccurrences[i].maxdegree = 0;
    fm::Bbrccandidatelegsoccurrences[i].selfjoin = 0;
    lastself[i] = NOTID;
    fm::Bbrccandidatelegsoccurrences[i].frequency = 0;
  }

  fm::Bbrccloselegsoccsused = false; // we are lazy with the initialization of close leg arrays, as we may not need them at all in
                             // many cases
  for ( BbrcOccurrenceId i = 0; i < legoccurrences.size (); i++ ) {
    BbrcLegOccurrence &legocc = legoccurrences[i];
    BbrcDatabaseTreePtr tree = fm::database->trees[legocc.tid];
    BbrcDatabaseTreeNode &node = tree->nodes[legocc.tonodeid];
    for ( int j = 0; j < node.edges.size (); j++ ) {
      if ( node.edges[j].tonode != legocc.fromnodeid ) {
	BbrcEdgeLabel edgelabel = node.edges[j].edgelabel;
        int number = nocycle ( tree, node, node.edges[j].tonode, i, &legoccurrencesdata );
        if ( number == 0 ) {
	  if ( edgelabel >= minlabel && edgelabel != neglect ) {
            vector<BbrcLegOccurrence> &candidatelegsoccs = fm::Bbrccandidatelegsoccurrences[edgelabel].elements;
            if ( candidatelegsoccs.empty () )
  	      fm::Bbrccandidatelegsoccurrences[edgelabel].frequency++;
	    else {
	      if ( candidatelegsoccs.back ().tid != legocc.tid )
  	        fm::Bbrccandidatelegsoccurrences[edgelabel].frequency++;
	      if ( candidatelegsoccs.back ().occurrenceid == i &&
                lastself[edgelabel] != (int) legocc.tid ) {
                lastself[edgelabel] = legocc.tid;
                fm::Bbrccandidatelegsoccurrences[edgelabel].selfjoin++;
              }
            }
            candidatelegsoccs.push_back ( BbrcLegOccurrence ( legocc.tid, i, node.edges[j].tonode, legocc.tonodeid ) );
	    Bbrcsetmax ( fm::Bbrccandidatelegsoccurrences[edgelabel].maxdegree, fm::database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
	  }
        }
        else if ( number - 1 != fm::graphstate->nodes.back().edges[0].tonode ) {
          candidateBbrcCloseBbrcLegsAllocate ( number, legoccurrencesdata.number + 1 );

          vector<CloseBbrcLegOccurrence> &candidatelegsoccs = fm::candidatecloselegsoccs[number][edgelabel].elements;
          if ( !candidatelegsoccs.size () || candidatelegsoccs.back ().tid != legocc.tid )
	    fm::candidatecloselegsoccs[number][edgelabel].frequency++;
          candidatelegsoccs.push_back ( CloseBbrcLegOccurrence ( legocc.tid, i ) );
          Bbrcsetmax ( fm::Bbrccandidatelegsoccurrences[edgelabel].maxdegree, fm::database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
        }
      }
    }
  }
}

/*

class Counter {
  public:
  int number;
  Counter () { number = 0; }
  ~Counter () { cout << number << endl; }
} counter, counter2;

void sanityCheck ( BbrcLegOccurrencesPtr legoccurrencesptr ) {
  if ( !legoccurrencesptr )
    return;
  counter2.number++;
  for ( BbrcOccurrenceId i = 0; i < legoccurrencesptr->elements.size (); i++ ) {
    for ( BbrcOccurrenceId t = i + 1; t < legoccurrencesptr->elements.size () &&
                                  legoccurrencesptr->elements[t].occurrenceid == legoccurrencesptr->elements[i].occurrenceid; t++ )
      if ( legoccurrencesptr->elements[t].tonodeid == legoccurrencesptr->elements[i].tonodeid )
        cerr << "MULTIPLE OCCURRENCES FOR THE SAME OCCURRENCE ID OF THE SAME NODE: "
             << legoccurrencesptr->elements[t].tonodeid << " IN " << legoccurrencesptr->elements[t].tid << endl;
    BbrcLegOccurrence &legocc = legoccurrencesptr->elements[i];
    BbrcOccurrenceId j = i, k;
    BbrcLegOccurrencesPtr walk1 = legoccurrencesptr, walk2;
    int val = 0;
    while ( walk1 ) {
      k = walk1->elements[j].occurrenceid;
      walk2 = walk1->parent;
      while ( walk2 ) {
        if ( walk2->elements[k].tid != walk1->elements[j].tid ) {
          cerr << "INSANITY: TIDs DO NOT MATCH: " << walk2->elements[k].tid << " AND " << walk1->elements[j].tid << endl;
        }
        if ( walk2->elements[k].tonodeid == walk1->elements[j].tonodeid ) {
          cerr << "INSANITY: SAME NODE TWICE: " << walk2->elements[k].tid << " (TID) " << walk1->elements[j].tonodeid << " (NODEID)" << endl;
        }
        k = walk2->elements[k].occurrenceid;
        walk2 = walk2->parent;
      }
      j = walk1->elements[j].occurrenceid;
      walk1 = walk1->parent;
      val++;
    }
    if ( val > counter.number )
      counter.number = val;
  }
}

*/
