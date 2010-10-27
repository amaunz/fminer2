// legoccurrence.cpp
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

#include "legoccurrence.h"
#include "closeleg.h"
#include "database.h"
#include "graphstate.h"

namespace fm {
    extern LastDatabase* last_database;
    extern LastGraphState* last_graphstate;
    extern LastLegOccurrences* last_legoccurrences;
    extern unsigned int last_minfreq;
    extern vector<LastLegOccurrences> last_Lastcandidatelegsoccurrences; 
    extern vector<vector< CloseLastLegOccurrences> > last_candidatecloselegsoccs;
    extern vector<bool> last_candidateLastcloselegsoccsused;
    extern bool last_Lastcloselegsoccsused;
}

void LastinitLastLegStatics () {
  fm::last_candidatecloselegsoccs.reserve ( 200 ); // should be larger than the largest structure that contains a cycle
  fm::last_Lastcandidatelegsoccurrences.resize ( fm::last_database->frequentLastEdgeLabelSize () );
}


/*
ostream &operator<< ( ostream &stream, LastLegOccurrence &occ ) {
  stream << "[" << occ.tid << "(" << fm::last_database->trees[occ.tid]->activity  << ")" << "," << occ.occurrenceid << "," << occ.tonodeid << "," << occ.fromnodeid << "]";
  return stream;
}
*/

ostream &operator<< ( ostream &stream, vector<LastLegOccurrence> &occs ) {
  LastTid lasttid = NOTID;
  LastFrequency frequency = 0;
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
LastLegOccurrencesPtr bbrc_join ( LastLegOccurrences &legoccsdata1, LastNodeId connectingnode, LastLegOccurrences &legoccsdata2 ) {
  if ( fm::last_graphstate->getNodeDegree ( connectingnode ) == fm::last_graphstate->getNodeMaxDegree ( connectingnode ) ) 
    return NULL;

  LastFrequency frequency = 0;
  LastTid lasttid = NOTID;
  vector<LastLegOccurrence> &legoccs1 = legoccsdata1.elements, &legoccs2 = legoccsdata2.elements;
  fm::last_legoccurrences->elements.resize ( 0 );
  fm::last_legoccurrences->maxdegree = 0;
  fm::last_legoccurrences->selfjoin = 0;
  //fm::last_legoccurrences->elements.reserve ( legoccs1.size () * 2 ); // increased memory usage, and speed!
  LastOccurrenceId j = 0, k = 0, l, m;
  unsigned int legoccs1size = legoccs1.size (), legoccs2size = legoccs2.size (); // this increases speed CONSIDERABLY!
  LastTid lastself = NOTID;

  do {
    while ( j < legoccs1size && legoccs1[j].occurrenceid < legoccs2[k].occurrenceid ) {
      j++;
    }
    if ( j < legoccs1size ) {
      LastLegOccurrence &jlegocc = legoccs1[j];
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
          for ( LastOccurrenceId m2 = m; m2 < j; m2++ ) {
            int d = 0;
            for ( LastOccurrenceId l2 = l; l2 < k; l2++ ) {
	      LastNodeId tonodeid = legoccs2[l2].tonodeid;
              if ( legoccs1[m2].tonodeid !=  tonodeid ) {
                fm::last_legoccurrences->elements.push_back ( LastLegOccurrence ( jlegocc.tid, m2, tonodeid, legoccs2[l2].fromnodeid ) );
                Lastsetmax ( fm::last_legoccurrences->maxdegree, fm::last_database->trees[jlegocc.tid]->nodes[tonodeid].edges.size () );
        		add = true;
        		d++;
              }
            }
	    if ( d > 1 && jlegocc.tid != lastself ) {
	      lastself = jlegocc.tid;
	      fm::last_legoccurrences->selfjoin++;
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

  if ( frequency >= fm::last_minfreq ) {
    fm::last_legoccurrences->parent = &legoccsdata1;
    fm::last_legoccurrences->number = legoccsdata1.number + 1;
    fm::last_legoccurrences->frequency = frequency;
    return fm::last_legoccurrences;
  }
  else
    return NULL;
}

LastLegOccurrencesPtr bbrc_join ( LastLegOccurrences &legoccsdata ) {
  if ( legoccsdata.selfjoin < fm::last_minfreq ) 
    return NULL;
  fm::last_legoccurrences->elements.resize ( 0 );
  vector<LastLegOccurrence> &legoccs = legoccsdata.elements;
  fm::last_legoccurrences->maxdegree = 0;
  fm::last_legoccurrences->selfjoin = 0;
  LastTid lastself = NOTID;

  LastOccurrenceId j = 0, k, l, m;
  do {
    k = j;
    LastLegOccurrence &legocc = legoccs[k];
    do {
      j++;
    }
    while ( j < legoccs.size () &&
            legoccs[j].occurrenceid == legocc.occurrenceid );
    for ( l = k; l < j; l++ )
      for ( m = k; m < j; m++ )
        if ( l != m ) {
          fm::last_legoccurrences->elements.push_back ( LastLegOccurrence ( legocc.tid, l, legoccs[m].tonodeid, legoccs[m].fromnodeid ) );
          Lastsetmax ( fm::last_legoccurrences->maxdegree, fm::last_database->trees[legocc.tid]->nodes[legoccs[m].tonodeid].edges.size () );
        }
    if ( ( j - k > 2 ) && legocc.tid != lastself ) {
      lastself = legocc.tid;
      fm::last_legoccurrences->selfjoin++;
    }
  }
  while ( j < legoccs.size () );

    // no need to check that we are frequent, we must be frequent
  fm::last_legoccurrences->parent = &legoccsdata;
  fm::last_legoccurrences->number = legoccsdata.number + 1;
  fm::last_legoccurrences->frequency = legoccsdata.selfjoin; 
    // we compute the self-join frequency exactly while building the
    // previous list. It is therefore not necessary to recompute it.
  return fm::last_legoccurrences;
}

inline int nocycle ( LastDatabaseTreePtr tree, LastDatabaseTreeNode &node, LastNodeId tonode, LastOccurrenceId occurrenceid, LastLegOccurrencesPtr legoccurrencesdataptr ) {
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

void candidateLastCloseLastLegsAllocate ( int number, int maxnumber ) {
  if ( !fm::last_Lastcloselegsoccsused ) {
    int oldsize = fm::last_candidatecloselegsoccs.size ();
    fm::last_candidatecloselegsoccs.resize ( maxnumber );
    for ( int k = oldsize; k < (int) fm::last_candidatecloselegsoccs.size (); k++ ) {
      fm::last_candidatecloselegsoccs[k].resize ( fm::last_database->frequentLastEdgeLabelSize () );
    }
    fm::last_candidateLastcloselegsoccsused.resize ( 0 );
    fm::last_candidateLastcloselegsoccsused.resize ( maxnumber, false );
    fm::last_Lastcloselegsoccsused = true;
  }
  if ( !fm::last_candidateLastcloselegsoccsused[number] ) {
    fm::last_candidateLastcloselegsoccsused[number] = true;
    vector<CloseLastLegOccurrences> &candidateedgelabeloccs = fm::last_candidatecloselegsoccs[number];
    for ( int k = 0; k < (int) candidateedgelabeloccs.size (); k++ ) {
      candidateedgelabeloccs[k].elements.resize ( 0 );
      candidateedgelabeloccs[k].frequency = 0;
    }
  }
}






void bbrc_extend ( LastLegOccurrences &legoccurrencesdata ) {
  // we're trying hard to avoid repeated destructor/constructor calls for complex types like vectors.
  // better reuse previously allocated memory, if possible!
  
  

  vector<LastLegOccurrence> &legoccurrences = legoccurrencesdata.elements;   ///////////////////////////////////////AM : BUG!!!



  LastTid lastself[fm::last_Lastcandidatelegsoccurrences.size ()];

  for ( int i = 0; i < (int) fm::last_Lastcandidatelegsoccurrences.size (); i++ ) {
    fm::last_Lastcandidatelegsoccurrences[i].elements.resize ( 0 );
    //fm::last_Lastcandidatelegsoccurrences[i].elements.reserve ( legoccurrences.size () ); // increases memory usage, but also speed!
    fm::last_Lastcandidatelegsoccurrences[i].parent = &legoccurrencesdata;
    fm::last_Lastcandidatelegsoccurrences[i].number = legoccurrencesdata.number + 1;
    fm::last_Lastcandidatelegsoccurrences[i].maxdegree = 0;
    fm::last_Lastcandidatelegsoccurrences[i].frequency = 0;
    fm::last_Lastcandidatelegsoccurrences[i].selfjoin = 0;
    lastself[i] = NOTID;
  }

  fm::last_Lastcloselegsoccsused = false; // we are lazy with the initialization of close leg arrays, as we may not need them at all in
                             // many cases

  for ( LastOccurrenceId i = 0; i < legoccurrences.size (); i++ ) {
    LastLegOccurrence &legocc = legoccurrences[i];
    LastDatabaseTreePtr tree = fm::last_database->trees[legocc.tid];
    LastDatabaseTreeNode &node = tree->nodes[legocc.tonodeid];
    for ( int j = 0; j < node.edges.size (); j++ ) {
      if ( node.edges[j].tonode != legocc.fromnodeid ) {
      	LastEdgeLabel edgelabel = node.edges[j].edgelabel;

        int number = nocycle ( tree, node, node.edges[j].tonode, i, &legoccurrencesdata );

        if ( number == 0 ) {
          vector<LastLegOccurrence> &candidatelegsoccs = fm::last_Lastcandidatelegsoccurrences[edgelabel].elements;
          if ( candidatelegsoccs.empty () )  fm::last_Lastcandidatelegsoccurrences[edgelabel].frequency++;
          else {

	            if ( candidatelegsoccs.back ().tid != legocc.tid )
        	        fm::last_Lastcandidatelegsoccurrences[edgelabel].frequency++;

	            if ( candidatelegsoccs.back ().occurrenceid == i &&
	                lastself[edgelabel] != legocc.tid ) {
                    lastself[edgelabel] = legocc.tid;
	                fm::last_Lastcandidatelegsoccurrences[edgelabel].selfjoin++;
	            }

          }
          candidatelegsoccs.push_back ( LastLegOccurrence ( legocc.tid, i, node.edges[j].tonode, legocc.tonodeid ) );
          Lastsetmax ( fm::last_Lastcandidatelegsoccurrences[edgelabel].maxdegree, fm::last_database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
        }

        else if ( number - 1 != fm::last_graphstate->nodes.back().edges[0].tonode ) {
            candidateLastCloseLastLegsAllocate ( number, legoccurrencesdata.number + 1 );
            vector<CloseLastLegOccurrence> &candidatelegsoccs = fm::last_candidatecloselegsoccs[number][edgelabel].elements;
            if ( !candidatelegsoccs.size () || candidatelegsoccs.back ().tid != legocc.tid )
	            fm::last_candidatecloselegsoccs[number][edgelabel].frequency++;
            candidatelegsoccs.push_back ( CloseLastLegOccurrence ( legocc.tid, i ) );
            Lastsetmax ( fm::last_Lastcandidatelegsoccurrences[edgelabel].maxdegree, fm::last_database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
        }

      }
    }
  }
}







void bbrc_extend ( LastLegOccurrences &legoccurrencesdata, LastEdgeLabel minlabel, LastEdgeLabel neglect ) {


  // we're trying hard to avoid repeated destructor/constructor calls for complex types like vectors.
  // better reuse previously allocated memory, if possible!




  vector<LastLegOccurrence> &legoccurrences = legoccurrencesdata.elements;  ///////////////////////////////////////AM : BUG!!!




  int lastself[fm::last_Lastcandidatelegsoccurrences.size ()];
  
  for ( int i = 0; i < (int) fm::last_Lastcandidatelegsoccurrences.size (); i++ ) {
    fm::last_Lastcandidatelegsoccurrences[i].elements.resize ( 0 );
    fm::last_Lastcandidatelegsoccurrences[i].parent = &legoccurrencesdata;
    fm::last_Lastcandidatelegsoccurrences[i].number = legoccurrencesdata.number + 1;
    fm::last_Lastcandidatelegsoccurrences[i].maxdegree = 0;
    fm::last_Lastcandidatelegsoccurrences[i].selfjoin = 0;
    lastself[i] = NOTID;
    fm::last_Lastcandidatelegsoccurrences[i].frequency = 0;
  }

  fm::last_Lastcloselegsoccsused = false; // we are lazy with the initialization of close leg arrays, as we may not need them at all in
                             // many cases
  for ( LastOccurrenceId i = 0; i < legoccurrences.size (); i++ ) {
    LastLegOccurrence &legocc = legoccurrences[i];
    LastDatabaseTreePtr tree = fm::last_database->trees[legocc.tid];
    LastDatabaseTreeNode &node = tree->nodes[legocc.tonodeid];
    for ( int j = 0; j < node.edges.size (); j++ ) {
      if ( node.edges[j].tonode != legocc.fromnodeid ) {
	LastEdgeLabel edgelabel = node.edges[j].edgelabel;
        int number = nocycle ( tree, node, node.edges[j].tonode, i, &legoccurrencesdata );
        if ( number == 0 ) {
	  if ( edgelabel >= minlabel && edgelabel != neglect ) {
            vector<LastLegOccurrence> &candidatelegsoccs = fm::last_Lastcandidatelegsoccurrences[edgelabel].elements;
            if ( candidatelegsoccs.empty () )
  	      fm::last_Lastcandidatelegsoccurrences[edgelabel].frequency++;
	    else {
	      if ( candidatelegsoccs.back ().tid != legocc.tid )
  	        fm::last_Lastcandidatelegsoccurrences[edgelabel].frequency++;
	      if ( candidatelegsoccs.back ().occurrenceid == i &&
                lastself[edgelabel] != (int) legocc.tid ) {
                lastself[edgelabel] = legocc.tid;
                fm::last_Lastcandidatelegsoccurrences[edgelabel].selfjoin++;
              }
            }
            candidatelegsoccs.push_back ( LastLegOccurrence ( legocc.tid, i, node.edges[j].tonode, legocc.tonodeid ) );
	    Lastsetmax ( fm::last_Lastcandidatelegsoccurrences[edgelabel].maxdegree, fm::last_database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
	  }
        }
        else if ( number - 1 != fm::last_graphstate->nodes.back().edges[0].tonode ) {
          candidateLastCloseLastLegsAllocate ( number, legoccurrencesdata.number + 1 );

          vector<CloseLastLegOccurrence> &candidatelegsoccs = fm::last_candidatecloselegsoccs[number][edgelabel].elements;
          if ( !candidatelegsoccs.size () || candidatelegsoccs.back ().tid != legocc.tid )
	    fm::last_candidatecloselegsoccs[number][edgelabel].frequency++;
          candidatelegsoccs.push_back ( CloseLastLegOccurrence ( legocc.tid, i ) );
          Lastsetmax ( fm::last_Lastcandidatelegsoccurrences[edgelabel].maxdegree, fm::last_database->trees[legocc.tid]->nodes[node.edges[j].tonode].edges.size () );
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

void sanityCheck ( LastLegOccurrencesPtr legoccurrencesptr ) {
  if ( !legoccurrencesptr )
    return;
  counter2.number++;
  for ( LastOccurrenceId i = 0; i < legoccurrencesptr->elements.size (); i++ ) {
    for ( LastOccurrenceId t = i + 1; t < legoccurrencesptr->elements.size () &&
                                  legoccurrencesptr->elements[t].occurrenceid == legoccurrencesptr->elements[i].occurrenceid; t++ )
      if ( legoccurrencesptr->elements[t].tonodeid == legoccurrencesptr->elements[i].tonodeid )
        cerr << "MULTIPLE OCCURRENCES FOR THE SAME OCCURRENCE ID OF THE SAME NODE: "
             << legoccurrencesptr->elements[t].tonodeid << " IN " << legoccurrencesptr->elements[t].tid << endl;
    LastLegOccurrence &legocc = legoccurrencesptr->elements[i];
    LastOccurrenceId j = i, k;
    LastLegOccurrencesPtr walk1 = legoccurrencesptr, walk2;
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
