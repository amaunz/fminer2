// path.cpp
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

#include <algorithm>
#include "patterntree.h"
#include "path.h"
#include "graphstate.h"
#include <iomanip>
#include "misc.h"

namespace fm {
    extern unsigned int bbrc_minfreq;
    extern bool bbrc_adjust_ub;
    extern bool bbrc_do_pruning;
    extern bool bbrc_do_backbone;
    extern bool bbrc_updated;
    extern int bbrc_type;
    extern bool bbrc_console_out;
    extern bool bbrc_refine_singles;
    extern bool bbrc_do_output;
    extern bool bbrc_bbrc_sep;
    extern bool bbrc_regression;

    extern BbrcDatabase* bbrc_database;
    extern ChisqBbrcConstraint* bbrc_chisq;
    extern KSBbrcConstraint* bbrc_ks;
    extern vector<string>* bbrc_result;
    extern BbrcStatistics* bbrc_statistics;
    extern BbrcGraphState* bbrc_graphstate;

    extern vector<BbrcLegOccurrences> bbrc_Bbrccandidatelegsoccurrences; 
}

// for every database node...
BbrcPath::BbrcPath ( BbrcNodeLabel startnodelabel ) {
  
    fm::bbrc_graphstate->insertStartNode ( startnodelabel );
    nodelabels.push_back ( startnodelabel );
    frontsymmetry = backsymmetry = totalsymmetry = 0;

    InputBbrcNodeLabel inl = fm::bbrc_database->nodelabels[startnodelabel].inputlabel;
    cerr << "Root: " << inl << endl;

    BbrcDatabaseBbrcNodeLabel &databasenodelabel = fm::bbrc_database->nodelabels[startnodelabel];

    // ...gather frequent edge labels
    vector<BbrcEdgeLabel> frequentedgelabels;
    for ( unsigned int i = 0; i < databasenodelabel.frequentedgelabels.size (); i++ )
        frequentedgelabels.push_back ( fm::bbrc_database->edgelabels[databasenodelabel.frequentedgelabels[i]].edgelabel );
                                                                                                    //  ^^^^^^^^^ is frequency rank!
    sort ( frequentedgelabels.begin (), frequentedgelabels.end () );                                // restores the rank order
    
    BbrcTid lastself[frequentedgelabels.size ()];
    vector<BbrcEdgeLabel> edgelabelorder ( fm::bbrc_database->edgelabelsindexes.size () );
    BbrcEdgeLabel j = 0;

    // FOR ALL EDGES...
    for ( unsigned int i = 0; i < frequentedgelabels.size (); i++ ) {
        edgelabelorder[frequentedgelabels[i]] = j;                      // For each rank, store it's position
        j++;
    
        // ...CREATE LEGS
        BbrcPathBbrcLegPtr leg = new BbrcPathBbrcLeg;
        legs.push_back ( leg );

        leg->tuple.depth = 0;                                           // TUPLE  DESCRIBES STRUCTURE...
        leg->tuple.connectingnode = 0;
        leg->tuple.edgelabel = frequentedgelabels[i]; 

        leg->occurrences.parent = &databasenodelabel.occurrences;       // ... OCCURRENCES DESCRIBES LOCATION IN TREE (1)
        leg->occurrences.number = 2;
        leg->occurrences.maxdegree = 0;
        leg->occurrences.selfjoin = 0;

        BbrcDatabaseBbrcEdgeLabel &databaseedgelabel = fm::bbrc_database->edgelabels[fm::bbrc_database->edgelabelsindexes[frequentedgelabels[i]]];
        leg->occurrences.frequency = databaseedgelabel.frequency;

        if ( databaseedgelabel.fromnodelabel == startnodelabel ) {
            leg->tuple.nodelabel = databaseedgelabel.tonodelabel;
        }

        else { 
            leg->tuple.nodelabel = databaseedgelabel.fromnodelabel;
        }

        lastself[i] = NOTID;

   }
    
    // ... OCCURRENCES DESCRIBES LOCATION IN TREE (2)
    for ( unsigned int i = 0; i < databasenodelabel.occurrences.elements.size (); i++ ) {
        BbrcDatabaseTree &tree = * (fm::bbrc_database->trees[databasenodelabel.occurrences.elements[i].tid]);
        BbrcDatabaseTreeNode &datanode = tree.nodes[databasenodelabel.occurrences.elements[i].tonodeid];
        for ( int j = 0; j < datanode.edges.size (); j++ ) {
            BbrcEdgeLabel edgelabel = edgelabelorder[datanode.edges[j].edgelabel];
            BbrcPathBbrcLeg &leg = * ( legs[edgelabel] );
            if ( !leg.occurrences.elements.empty () &&
                  leg.occurrences.elements.back ().occurrenceid == i &&
                  lastself[edgelabel] != tree.tid ) {
                //leg.occurrences.selfjoin++;
                leg.occurrences.selfjoin += 
                fm::bbrc_database->trees[tree.tid]->weight;
                lastself[edgelabel] = tree.tid;
            }
            vector_push_back ( BbrcLegOccurrence, leg.occurrences.elements, legoccurrence );
            legoccurrence.tid = tree.tid;
            legoccurrence.occurrenceid = i;
            legoccurrence.tonodeid = datanode.edges[j].tonode;
            legoccurrence.fromnodeid = databasenodelabel.occurrences.elements[i].tonodeid;
        }
    }
  
}

BbrcPath::BbrcPath ( BbrcPath &parentpath, unsigned int legindex ) {
  BbrcPathBbrcLeg &leg = (*parentpath.legs[legindex]);
  int positionshift;
  
  // fill in normalisation information, it seems a lot of code, but in fact it's just a walk through the edge/nodelabels arrays.

  nodelabels.resize ( parentpath.nodelabels.size () + 1 );
  edgelabels.resize ( parentpath.edgelabels.size () + 1 );

  BbrcaddCloseExtensions ( closelegs, parentpath.closelegs, leg.occurrences );

  if ( parentpath.nodelabels.size () == 1 ) {
    totalsymmetry = parentpath.nodelabels[0] - leg.tuple.nodelabel;
    frontsymmetry = backsymmetry = 0;
    nodelabels[1] = leg.tuple.nodelabel;
    edgelabels[0] = leg.tuple.edgelabel;
    nodelabels[0] = parentpath.nodelabels[0];
    positionshift = 0;
  }
  else if ( leg.tuple.depth == 0 ) {
    positionshift = 1;
    nodelabels[0] = leg.tuple.nodelabel;
    edgelabels[0] = leg.tuple.edgelabel;

    backsymmetry = parentpath.totalsymmetry;
    frontsymmetry = leg.tuple.nodelabel - parentpath.nodelabels[parentpath.nodelabels.size () - 2];
    totalsymmetry = leg.tuple.nodelabel - parentpath.nodelabels.back ();
    if ( !totalsymmetry )
      totalsymmetry = leg.tuple.edgelabel - parentpath.edgelabels.back ();

    unsigned int i = 0;
    // we can prepend only before strings of length 2
    if ( parentpath.nodelabels.size () > 2 ) {
      if ( !frontsymmetry )
        frontsymmetry = leg.tuple.edgelabel - parentpath.edgelabels[parentpath.nodelabels.size () - 3];

      while ( !frontsymmetry && i < parentpath.edgelabels.size () / 2 ) {
        nodelabels[i + 1] = parentpath.nodelabels[i];
        edgelabels[i + 1] = parentpath.edgelabels[i];

        frontsymmetry = parentpath.nodelabels[i] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 3];
        if ( !frontsymmetry && parentpath.nodelabels.size () > 3 )
          frontsymmetry = parentpath.edgelabels[i] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 4];

	    if ( !totalsymmetry ) {
	        totalsymmetry = parentpath.nodelabels[i] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 2];
	    if ( !totalsymmetry )
	        totalsymmetry = parentpath.edgelabels[i] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 3];
	    }

        i++;
      }
    }
    for ( ; !totalsymmetry && i < parentpath.edgelabels.size () / 2; i++ ) {
      nodelabels[i + 1] = parentpath.nodelabels[i];
      edgelabels[i + 1] = parentpath.edgelabels[i];

      totalsymmetry = parentpath.nodelabels[i] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 2];
      if ( !totalsymmetry && parentpath.nodelabels.size () > 3 )
        totalsymmetry = parentpath.edgelabels[i] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 3];
    }
    for ( ;i < parentpath.edgelabels.size (); i++ ) {
      nodelabels[i + 1] = parentpath.nodelabels[i];
      edgelabels[i + 1] = parentpath.edgelabels[i];
    }

    nodelabels[i + 1] = parentpath.nodelabels[i];



    // build OccurrenceLists
    bbrc_extend ( leg.occurrences );
    for (unsigned int i = 0; i < fm::bbrc_Bbrccandidatelegsoccurrences.size (); i++ ) {
      if ( fm::bbrc_Bbrccandidatelegsoccurrences[i].frequency >= fm::bbrc_minfreq ) {
        BbrcPathBbrcLegPtr leg2 = new BbrcPathBbrcLeg;
        legs.push_back ( leg2 );
        leg2->tuple.edgelabel = i;
    	leg2->tuple.connectingnode = fm::bbrc_graphstate->lastNode ();
        BbrcDatabaseBbrcEdgeLabel &databaseedgelabel = fm::bbrc_database->edgelabels[fm::bbrc_database->edgelabelsindexes[i]];
        if ( databaseedgelabel.fromnodelabel == leg.tuple.nodelabel )
          leg2->tuple.nodelabel = databaseedgelabel.tonodelabel;
        else
          leg2->tuple.nodelabel = databaseedgelabel.fromnodelabel;
        leg2->tuple.depth = 0;
        store ( leg2->occurrences, fm::bbrc_Bbrccandidatelegsoccurrences[i] ); // avoid copying
      }
    }



  }


  else  {
    positionshift = 0;

    frontsymmetry = parentpath.totalsymmetry;
    backsymmetry = parentpath.nodelabels[1] - leg.tuple.nodelabel;
    totalsymmetry = parentpath.nodelabels[0] - leg.tuple.nodelabel;
    if ( !totalsymmetry )
      totalsymmetry = parentpath.edgelabels[0] - leg.tuple.edgelabel;
    unsigned int i = 0;
    if ( parentpath.nodelabels.size () > 2 ) {
      if ( !backsymmetry )
        backsymmetry = parentpath.edgelabels[1] - leg.tuple.edgelabel;

      while ( !backsymmetry && i < parentpath.edgelabels.size () / 2 ) {
        nodelabels[i] = parentpath.nodelabels[i];
	edgelabels[i] = parentpath.edgelabels[i];

	backsymmetry = parentpath.nodelabels[i + 2] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 1];
	if ( !backsymmetry && parentpath.nodelabels.size () > 3 )
	  backsymmetry = parentpath.edgelabels[i + 2] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 2];

	if ( !totalsymmetry ) {
	  totalsymmetry = parentpath.nodelabels[i + 1] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 1];
	  if ( !totalsymmetry && parentpath.nodelabels.size () > 3 )
	    totalsymmetry = parentpath.edgelabels[i + 1] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 2];
	}
	i++;
      }
    }
    for ( ; !totalsymmetry && i < parentpath.edgelabels.size () / 2; i++ ) {
      nodelabels[i] = parentpath.nodelabels[i];
      edgelabels[i] = parentpath.edgelabels[i];
      totalsymmetry = parentpath.nodelabels[i + 1] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 1];
      if ( !totalsymmetry && i < parentpath.edgelabels.size () - 1 )
        totalsymmetry = parentpath.edgelabels[i + 1] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 2];
    }
    for ( ; i < parentpath.edgelabels.size (); i++ ) {
      nodelabels[i] = parentpath.nodelabels[i];
      edgelabels[i] = parentpath.edgelabels[i];
    }

    nodelabels[i] = parentpath.nodelabels[i];
    edgelabels[i] = leg.tuple.edgelabel;
    nodelabels[i+1] = leg.tuple.nodelabel;
  }

  unsigned int i = 0;
  BbrcLegOccurrencesPtr legoccurrencesptr;
  for ( ; i < legindex; i++ ) {
    BbrcPathBbrcLeg &leg2 = (*parentpath.legs[i]);

    if ( (legoccurrencesptr = bbrc_join ( leg.occurrences, leg2.tuple.connectingnode, leg2.occurrences )) ) { // JOIN OCCURRENCES
      BbrcPathBbrcLegPtr leg3 = new BbrcPathBbrcLeg;
      legs.push_back ( leg3 );
      leg3->tuple.connectingnode = leg2.tuple.connectingnode;
      leg3->tuple.edgelabel = leg2.tuple.edgelabel;
      leg3->tuple.nodelabel = leg2.tuple.nodelabel;
      leg3->tuple.depth = leg2.tuple.depth + positionshift;
      store ( leg3->occurrences, *legoccurrencesptr );
    }
  }

  if ( (legoccurrencesptr = bbrc_join ( leg.occurrences )) ) {
    BbrcPathBbrcLegPtr leg3 = new BbrcPathBbrcLeg;
    legs.push_back ( leg3 );
    leg3->tuple.connectingnode = leg.tuple.connectingnode;
    leg3->tuple.edgelabel = leg.tuple.edgelabel;
    leg3->tuple.nodelabel = leg.tuple.nodelabel;
    leg3->tuple.depth = leg.tuple.depth + positionshift;
    store ( leg3->occurrences, *legoccurrencesptr );
  }

  for ( i++; i < parentpath.legs.size (); i++ ) {
    BbrcPathBbrcLeg &leg2 = (*parentpath.legs[i]);
    if ( (legoccurrencesptr = bbrc_join ( leg.occurrences, leg2.tuple.connectingnode, leg2.occurrences )) ) {
      BbrcPathBbrcLegPtr leg3 = new BbrcPathBbrcLeg;
      legs.push_back ( leg3 );
      leg3->tuple.connectingnode = leg2.tuple.connectingnode;
      leg3->tuple.edgelabel = leg2.tuple.edgelabel;
      leg3->tuple.nodelabel = leg2.tuple.nodelabel;
      leg3->tuple.depth = leg2.tuple.depth + positionshift;
      store ( leg3->occurrences, *legoccurrencesptr );
    }
  }

  if ( positionshift ) {
    BbrcaddCloseExtensions ( closelegs, leg.occurrences.number ); // stored separately
    return;
  }

  bbrc_extend ( leg.occurrences );
  for ( unsigned int i = 0; i < fm::bbrc_Bbrccandidatelegsoccurrences.size (); i++ ) {
    if ( fm::bbrc_Bbrccandidatelegsoccurrences[i].frequency >= fm::bbrc_minfreq ) {
      BbrcPathBbrcLegPtr leg2 = new BbrcPathBbrcLeg;
      legs.push_back ( leg2 );
      leg2->tuple.edgelabel = i;
      leg2->tuple.connectingnode = fm::bbrc_graphstate->lastNode ();
      BbrcDatabaseBbrcEdgeLabel &databaseedgelabel = fm::bbrc_database->edgelabels[fm::bbrc_database->edgelabelsindexes[i]];
      if ( databaseedgelabel.fromnodelabel == leg.tuple.nodelabel )
        leg2->tuple.nodelabel = databaseedgelabel.tonodelabel;
      else
        leg2->tuple.nodelabel = databaseedgelabel.fromnodelabel;
      leg2->tuple.depth = leg.tuple.depth + 1;
      store ( leg2->occurrences, fm::bbrc_Bbrccandidatelegsoccurrences[i] ); // avoid copying
    }
  }

  BbrcaddCloseExtensions ( closelegs, leg.occurrences.number );
}

BbrcPath::~BbrcPath () {
  for ( unsigned int i = 0; i < legs.size (); i++ )
    delete legs[i];
  for ( unsigned int i = 0; i < closelegs.size (); i++ )
    delete closelegs[i];
}

// ADDED
bool BbrcPath::is_normal ( BbrcEdgeLabel edgelabel ) {
  // symplistic quadratic algorithm
  int nodelabelssize = nodelabels.size (), step, add, start;
  
  edgelabels.push_back ( edgelabel );
    
  // if we would program it better, we would use the 'totalsymmetry' variable here;
  // however, to be quick & easy, we used a different coding here...
  int t = nodelabelssize - 1, r = 0;
  int symmetry;
  do {
    symmetry = nodelabels[t] - nodelabels[r];
    int nt = ( t + nodelabelssize - 1 ) % nodelabelssize;
    if ( !symmetry ) 
      symmetry = edgelabels[nt] - edgelabels[r];
    r = ( r + 1 ) % nodelabelssize;
    t = nt;
  }
  while ( symmetry == 0 && t != nodelabelssize - 1 );
  
  if ( symmetry < 0 ) {
    step = -1 + nodelabelssize ;
    add = -1 + nodelabelssize ;
    start = nodelabelssize - 1;
  }
  else {
    step = 1 + nodelabelssize;
    add = nodelabelssize ;
    start = 0;
  }
  for ( int i = 0; i < nodelabelssize; i++ ) {
    // starting positions for the new path
    int k = start, l = i, p;
    do {
      if ( nodelabels[l] < nodelabels[k] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( nodelabels[l] > nodelabels[k] )
        break;
      p = ( k + add ) % nodelabelssize;
      l = ( l + nodelabelssize - 1 ) % nodelabelssize;
      if ( edgelabels[l] < edgelabels[p] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( edgelabels[l] > edgelabels[p] ) 
        break;
      k = ( k + step ) % nodelabelssize;
    }
    while ( k != start );
    
    k = start, l = i;
    do {
      if ( nodelabels[l] < nodelabels[k] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( nodelabels[l] > nodelabels[k] ) 
        break;
      p = ( k + add ) % nodelabelssize;
      if ( edgelabels[l] < edgelabels[p] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( edgelabels[l] > edgelabels[p] ) 
        break;
      l = ( l + 1 ) % nodelabelssize;
      k = ( k + step ) % nodelabelssize;
    }
    while ( k != start );
    
  }
  edgelabels.pop_back ();
  return true;
}









void BbrcPath::expand2 (pair<float,string> max) {

  fm::bbrc_statistics->patternsize++;
  if ( (unsigned) fm::bbrc_statistics->patternsize > fm::bbrc_statistics->frequenttreenumbers.size () ) {
    fm::bbrc_statistics->frequenttreenumbers.push_back ( 0 );
    fm::bbrc_statistics->frequentpathnumbers.push_back ( 0 );
    fm::bbrc_statistics->frequentgraphnumbers.push_back ( 0 );
  }
  ++fm::bbrc_statistics->frequentpathnumbers[fm::bbrc_statistics->patternsize-1];
  
  if ( fm::bbrc_statistics->patternsize == ((1<<(sizeof(BbrcNodeId)*8))-1) ) {
    fm::bbrc_statistics->patternsize--;
    return;
  }

  vector<unsigned int> forwpathlegs; forwpathlegs.clear();
  vector<unsigned int> backwpathlegs; backwpathlegs.clear();
  vector<unsigned int> pathlegs; pathlegs.clear();

  for ( unsigned int i = 0; i < legs.size (); i++ ) {

    BbrcPathBbrcTuple &tuple = legs[i]->tuple;
    if ( tuple.depth == nodelabels.size () - 1 ) {

      if ( tuple.nodelabel > nodelabels[0] ||
           ( tuple.nodelabel == nodelabels[0] &&
             ( tuple.edgelabel > edgelabels[0] ||
               ( tuple.edgelabel == edgelabels[0] && backsymmetry <= 0 )
             )
           ) ) {
        forwpathlegs.push_back(i);
        pathlegs.push_back(i);

      }
    }
  }

  for ( unsigned int i = 0; i < legs.size (); i++ ) {
    BbrcPathBbrcTuple &tuple = legs[i]->tuple;
    if ( tuple.depth != nodelabels.size () - 1 ) {


      if ( legs[i]->tuple.depth == 0 ) {
        if ( totalsymmetry &&
             ( tuple.nodelabel > nodelabels.back () ||
             ( tuple.nodelabel == nodelabels.back () &&
               ( tuple.edgelabel > edgelabels.back () ||
                 ( tuple.edgelabel == edgelabels.back () && frontsymmetry >= 0 )
               )
             ) ) ) {
            backwpathlegs.push_back(i);
            pathlegs.push_back(i);
        }
      }
    }
  }
 

  // we have reached a leaf
  if (fm::bbrc_do_backbone && (pathlegs.size()==0)) { 
     if (fm::bbrc_updated) { 
        if (fm::bbrc_do_output) {
            if (!fm::bbrc_console_out) { 
                (*fm::bbrc_result) << max.second; 
            }
            else cout << max.second;
        }
        fm::bbrc_updated = false;
     }
  }

  
  
  
  // Grow BbrcPath forw
  for (unsigned int j=0; j<forwpathlegs.size() ; j++ ) {
    unsigned int index = forwpathlegs[j];


    // Calculate chisq
     if (fm::bbrc_chisq->active) { 
        if (!fm::bbrc_regression) fm::bbrc_chisq->Calc(legs[index]->occurrences.elements);
        else fm::bbrc_ks->Calc(legs[index]->occurrences.elements);
      }
          
    // GRAPHSTATE AND OUTPUT
    fm::bbrc_graphstate->insertNode ( legs[index]->tuple.connectingnode, legs[index]->tuple.edgelabel, legs[index]->occurrences.maxdegree );


    // immediate output
    if (fm::bbrc_do_output && !fm::bbrc_do_backbone) {
        if (!fm::bbrc_console_out) (*fm::bbrc_result) << fm::bbrc_graphstate->to_s(legs[index]->occurrences.frequency);
        else fm::bbrc_graphstate->print(legs[index]->occurrences.frequency);
    }


    // RECURSE
    float cmax = maxi ( maxi ( fm::bbrc_chisq->sig, max.first ), fm::bbrc_chisq->p );
    if ( (
             !fm::bbrc_do_pruning || 
             (
               (  !fm::bbrc_adjust_ub && (fm::bbrc_chisq->u >= fm::bbrc_chisq->sig) ) || 
               (   fm::bbrc_adjust_ub && (fm::bbrc_chisq->u >= cmax) )
             )
         ) 
            &&
         (
            fm::bbrc_refine_singles || (legs[index]->occurrences.frequency>1)
         )
      ){   // UB-PRUNING

      BbrcPath path ( *this, index );
      if (!fm::bbrc_regression) {
          if (max.first<fm::bbrc_chisq->p) { fm::bbrc_updated = true; path.expand2 ( pair<float, string>(fm::bbrc_chisq->p, fm::bbrc_graphstate->to_s(legs[index]->occurrences.frequency))); }
          else path.expand2 (max);
      }
      else {
          if (max.first<fm::bbrc_ks->p) { fm::bbrc_updated = true; path.expand2 ( pair<float, string>(fm::bbrc_ks->p, fm::bbrc_graphstate->to_s(legs[index]->occurrences.frequency))); }
          else path.expand2 (max);
      }
    }
    else {
        if (fm::bbrc_do_backbone && fm::bbrc_updated) {  // FREE STRUCTURES: search was pruned
            if (fm::bbrc_do_output) {
                if (!fm::bbrc_console_out) (*fm::bbrc_result) << max.second;
                else cout << max.second;
            }
            fm::bbrc_updated=false;
        }
    }

    fm::bbrc_graphstate->deleteNode ();
  }



  // Grow BbrcPath backw
  for (unsigned int j=0; j<backwpathlegs.size() ; j++ ) {
    unsigned int index = backwpathlegs[j];
    

    // Calculate chisq
    if (fm::bbrc_chisq->active) { 
        if (!fm::bbrc_regression) fm::bbrc_chisq->Calc(legs[index]->occurrences.elements);
        else fm::bbrc_ks->Calc(legs[index]->occurrences.elements);
    }


    // GRAPHSTATE AND OUTPUT
    fm::bbrc_graphstate->insertNode ( legs[index]->tuple.connectingnode, legs[index]->tuple.edgelabel, legs[index]->occurrences.maxdegree );

    // immediate output
    if (fm::bbrc_do_output && !fm::bbrc_do_backbone) {
        if (!fm::bbrc_console_out) (*fm::bbrc_result) << fm::bbrc_graphstate->to_s(legs[index]->occurrences.frequency);
        else fm::bbrc_graphstate->print(legs[index]->occurrences.frequency);
    }

    // RECURSE
    float cmax = maxi ( maxi ( fm::bbrc_chisq->sig, max.first ), fm::bbrc_chisq->p );
    if ( ( 
             !fm::bbrc_do_pruning || 
             (
               (  !fm::bbrc_adjust_ub && (fm::bbrc_chisq->u >= fm::bbrc_chisq->sig) ) || 
               (   fm::bbrc_adjust_ub && (fm::bbrc_chisq->u >= cmax) )
             )
         ) 
            &&
         (
            fm::bbrc_refine_singles || (legs[index]->occurrences.frequency>1)
         )
     ){   // UB-PRUNING

      BbrcPath path ( *this, index );
      if (!fm::bbrc_regression) {
          if (max.first<fm::bbrc_chisq->p) { fm::bbrc_updated = true; path.expand2 ( pair<float, string>(fm::bbrc_chisq->p, fm::bbrc_graphstate->to_s(legs[index]->occurrences.frequency))); }
          else path.expand2 (max);
      }
      else {
          if (max.first<fm::bbrc_ks->p) { fm::bbrc_updated = true; path.expand2 ( pair<float, string>(fm::bbrc_ks->p, fm::bbrc_graphstate->to_s(legs[index]->occurrences.frequency))); }
          else path.expand2 (max);
      }
    }
    else {
        if (fm::bbrc_do_backbone && fm::bbrc_updated) { // FREE STRUCTURES: search was pruned
            if (fm::bbrc_do_output) {
                if (!fm::bbrc_console_out) (*fm::bbrc_result) << max.second;
                else if (fm::bbrc_do_output) cout << max.second;
            }
            fm::bbrc_updated=false;
        }
    }

    fm::bbrc_graphstate->deleteNode ();

  }



  bool uptmp = fm::bbrc_updated;

  if (fm::bbrc_bbrc_sep && !fm::bbrc_do_backbone && legs.size() > 0) {
      if (fm::bbrc_do_output && !fm::bbrc_console_out && fm::bbrc_result->size() && (fm::bbrc_result->back()!=fm::bbrc_graphstate->sep())) (*fm::bbrc_result) << fm::bbrc_graphstate->sep();
      //else cout << fm::bbrc_graphstate->sep() << endl;
  }

  for ( unsigned int i = 0; i < legs.size (); i++ ) {
    BbrcPathBbrcTuple &tuple = legs[i]->tuple;
    if ( tuple.depth != nodelabels.size () - 1 ) {


      // PHASE 2: GROW TREE
      if ( legs[i]->tuple.depth != 0 ) {
        if ( ( totalsymmetry || legs[i]->tuple.depth <= edgelabels.size () / 2 ) &&
 	      ( legs[i]->tuple.depth != 1 || legs[i]->tuple.edgelabel >= edgelabels[0] ) &&
	      ( legs[i]->tuple.depth != nodelabels.size () - 2 || legs[i]->tuple.edgelabel >= edgelabels.back () ) &&
	        fm::bbrc_type > 1 ) {
          // Calculate chisq
          if (fm::bbrc_chisq->active) { 
            if (!fm::bbrc_regression) fm::bbrc_chisq->Calc(legs[i]->occurrences.elements);
            else fm::bbrc_ks->Calc(legs[i]->occurrences.elements);
          }


          // GRAPHSTATE
          fm::bbrc_graphstate->insertNode ( legs[i]->tuple.connectingnode, legs[i]->tuple.edgelabel, legs[i]->occurrences.maxdegree );

          // immediate output
          if (fm::bbrc_do_output && !fm::bbrc_do_backbone) {
             if (!fm::bbrc_console_out) (*fm::bbrc_result) << fm::bbrc_graphstate->to_s(legs[i]->occurrences.frequency);
             else fm::bbrc_graphstate->print(legs[i]->occurrences.frequency);
          }

          // RECURSE
          float cmax = maxi ( maxi ( fm::bbrc_chisq->sig, max.first ), fm::bbrc_chisq->p );

          if ( ( !fm::bbrc_do_pruning || 
               (
                 (  !fm::bbrc_adjust_ub && (fm::bbrc_chisq->u >= fm::bbrc_chisq->sig) ) || 
                 (   fm::bbrc_adjust_ub && (fm::bbrc_chisq->u >= cmax) )
               )
             ) &&
             (
                fm::bbrc_refine_singles || (legs[i]->occurrences.frequency>1)
             )
    
          ){   // UB-PRUNING

            BbrcPatternTree tree ( *this, i );

            if (!fm::bbrc_regression) {
                if (max.first<fm::bbrc_chisq->p) { fm::bbrc_updated = true; tree.expand ( pair<float, string>(fm::bbrc_chisq->p, fm::bbrc_graphstate->to_s(legs[i]->occurrences.frequency))); }
                else tree.expand (max);
            }
            else {
                if (max.first<fm::bbrc_ks->p) { fm::bbrc_updated = true; tree.expand ( pair<float, string>(fm::bbrc_ks->p, fm::bbrc_graphstate->to_s(legs[i]->occurrences.frequency))); }
                else tree.expand (max);
            }

          }

          else {
            if (fm::bbrc_do_backbone && fm::bbrc_updated) { 
              if (fm::bbrc_do_output) {
                  if (!fm::bbrc_console_out) (*fm::bbrc_result) << max.second;
                  else  cout << max.second;
              }
              fm::bbrc_updated=false;
            }
          }

	      fm::bbrc_graphstate->deleteNode ();
        }
      }
    }
  }


  fm::bbrc_updated=uptmp;
    
  fm::bbrc_statistics->patternsize--;


//  cerr << "backtracking p" << endl;
}








void BbrcPath::expand () {

  for ( unsigned int i = 0; i < legs.size (); i++ ) {
    BbrcPathBbrcTuple &tuple = legs[i]->tuple;
    if ( tuple.nodelabel >= nodelabels[0] ) {
        
      if (fm::bbrc_chisq->active) { 
        if (!fm::bbrc_regression) fm::bbrc_chisq->Calc(legs[i]->occurrences.elements);
        else fm::bbrc_ks->Calc(legs[i]->occurrences.elements);
      }

      // GRAPHSTATE AND OUTPUT
      fm::bbrc_graphstate->insertNode ( tuple.connectingnode, tuple.edgelabel, legs[i]->occurrences.maxdegree );
      if (fm::bbrc_do_output && !fm::bbrc_do_backbone && legs[i]->occurrences.frequency>=fm::bbrc_minfreq) { 
          if (!fm::bbrc_console_out) (*fm::bbrc_result) << fm::bbrc_graphstate->to_s(legs[i]->occurrences.frequency);
          else fm::bbrc_graphstate->print(legs[i]->occurrences.frequency);
      }

      // RECURSE
      BbrcPath path (*this, i);
      fm::bbrc_updated = true;
      path.expand2 (pair<float, string>(fm::bbrc_chisq->p, fm::bbrc_graphstate->to_s(legs[i]->occurrences.frequency)));
      fm::bbrc_graphstate->deleteNode ();

    }
  }
  fm::bbrc_graphstate->deleteStartNode ();


//  cerr << "backtracking p" << endl;
}






ostream &operator<< ( ostream &stream, BbrcPath &path ) {
  stream << /* database->nodelabels[ */ (int) path.nodelabels[0] /* ].inputlabel; */ << " ";
  for ( unsigned int i = 0; i < path.edgelabels.size (); i++ ) {
    //stream << (char) ( path.edgelabels[i] + 'A' ) << path.nodelabels[i+1];
    stream << /*database->edgelabels[database->edgelabelsindexes[*/ (int) path.edgelabels[i] /*]].inputedgelabel */ << " " <<  /* database->nodelabels[ */ (int) path.nodelabels[i+1] /* ].inputlabel */ << " ";
  }
  return stream;
}
