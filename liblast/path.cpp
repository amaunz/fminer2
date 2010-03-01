// path.cpp
// © 2008 by Andreas Maunz, andreas@maunz.de, jul 2008
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.

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

#include <algorithm>
#include "patterntree.h"
#include "path.h"
#include "graphstate.h"
#include <iomanip>
#include "misc.h"

namespace fm {
    extern unsigned int minfreq;
    extern bool do_pruning;
    extern bool updated;
    extern int type;
    extern bool console_out;
    extern bool refine_singles;
    extern bool do_output;
    extern bool bbrc_sep;
    extern bool most_specific_trees_only;
    extern bool gsp_out;
    extern bool die;
    extern bool do_last;

    extern Database* database;
    extern ChisqConstraint* chisq;
    extern vector<string>* result;
    extern Statistics* statistics;
    extern GraphState* graphstate;

    extern vector<LegOccurrences> candidatelegsoccurrences; 
}

// for every database node...
Path::Path ( NodeLabel startnodelabel ) {
  
    fm::graphstate->insertStartNode ( startnodelabel );
    nodelabels.push_back ( startnodelabel );
    frontsymmetry = backsymmetry = totalsymmetry = 0;

    InputNodeLabel inl = fm::database->nodelabels[startnodelabel].inputlabel;
    cerr << "Root: " << inl << endl;

    DatabaseNodeLabel &databasenodelabel = fm::database->nodelabels[startnodelabel];

    // ...gather frequent edge labels
    vector<EdgeLabel> frequentedgelabels;
    for ( unsigned int i = 0; i < databasenodelabel.frequentedgelabels.size (); i++ )
        frequentedgelabels.push_back ( fm::database->edgelabels[databasenodelabel.frequentedgelabels[i]].edgelabel );
                                                                                                    //  ^^^^^^^^^ is frequency rank!
    sort ( frequentedgelabels.begin (), frequentedgelabels.end () );                                // restores the rank order
    
    Tid lastself[frequentedgelabels.size ()];
    vector<EdgeLabel> edgelabelorder ( fm::database->edgelabelsindexes.size () );
    EdgeLabel j = 0;

    // FOR ALL EDGES...
    for ( unsigned int i = 0; i < frequentedgelabels.size (); i++ ) {
        edgelabelorder[frequentedgelabels[i]] = j;                      // For each rank, store it's position
        j++;
    
        // ...CREATE LEGS
        PathLegPtr leg = new PathLeg;
        legs.push_back ( leg );

        leg->tuple.depth = 0;                                           // TUPLE  DESCRIBES STRUCTURE...
        leg->tuple.connectingnode = 0;
        leg->tuple.edgelabel = frequentedgelabels[i]; 

        leg->occurrences.parent = &databasenodelabel.occurrences;       // ... OCCURRENCES DESCRIBES LOCATION IN TREE (1)
        leg->occurrences.number = 2;
        leg->occurrences.maxdegree = 0;
        leg->occurrences.selfjoin = 0;

        DatabaseEdgeLabel &databaseedgelabel = fm::database->edgelabels[fm::database->edgelabelsindexes[frequentedgelabels[i]]];
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
        DatabaseTree &tree = * (fm::database->trees[databasenodelabel.occurrences.elements[i].tid]);
        DatabaseTreeNode &datanode = tree.nodes[databasenodelabel.occurrences.elements[i].tonodeid];
        for ( int j = 0; j < datanode.edges.size (); j++ ) {
            EdgeLabel edgelabel = edgelabelorder[datanode.edges[j].edgelabel];
            PathLeg &leg = * ( legs[edgelabel] );
            if ( !leg.occurrences.elements.empty () &&
                  leg.occurrences.elements.back ().occurrenceid == i &&
                  lastself[edgelabel] != tree.tid ) {
                leg.occurrences.selfjoin++;
                lastself[edgelabel] = tree.tid;
            }
            vector_push_back ( LegOccurrence, leg.occurrences.elements, legoccurrence );
            legoccurrence.tid = tree.tid;
            legoccurrence.occurrenceid = i;
            legoccurrence.tonodeid = datanode.edges[j].tonode;
            legoccurrence.fromnodeid = databasenodelabel.occurrences.elements[i].tonodeid;
        }
    }
  
}

Path::Path ( Path &parentpath, unsigned int legindex ) {
  PathLeg &leg = (*parentpath.legs[legindex]);
  int positionshift;
  
  // fill in normalisation information, it seems a lot of code, but in fact it's just a walk through the edge/nodelabels arrays.

  nodelabels.resize ( parentpath.nodelabels.size () + 1 );
  edgelabels.resize ( parentpath.edgelabels.size () + 1 );

  addCloseExtensions ( closelegs, parentpath.closelegs, leg.occurrences );

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
    extend ( leg.occurrences );
    for (unsigned int i = 0; i < fm::candidatelegsoccurrences.size (); i++ ) {
      if ( fm::candidatelegsoccurrences[i].frequency >= fm::minfreq ) {
        PathLegPtr leg2 = new PathLeg;
        legs.push_back ( leg2 );
        leg2->tuple.edgelabel = i;
    	leg2->tuple.connectingnode = fm::graphstate->lastNode ();
        DatabaseEdgeLabel &databaseedgelabel = fm::database->edgelabels[fm::database->edgelabelsindexes[i]];
        if ( databaseedgelabel.fromnodelabel == leg.tuple.nodelabel )
          leg2->tuple.nodelabel = databaseedgelabel.tonodelabel;
        else
          leg2->tuple.nodelabel = databaseedgelabel.fromnodelabel;
        leg2->tuple.depth = 0;
        store ( leg2->occurrences, fm::candidatelegsoccurrences[i] ); // avoid copying
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
  LegOccurrencesPtr legoccurrencesptr;
  for ( ; i < legindex; i++ ) {
    PathLeg &leg2 = (*parentpath.legs[i]);

    if ( (legoccurrencesptr = join ( leg.occurrences, leg2.tuple.connectingnode, leg2.occurrences )) ) { // JOIN OCCURRENCES
      PathLegPtr leg3 = new PathLeg;
      legs.push_back ( leg3 );
      leg3->tuple.connectingnode = leg2.tuple.connectingnode;
      leg3->tuple.edgelabel = leg2.tuple.edgelabel;
      leg3->tuple.nodelabel = leg2.tuple.nodelabel;
      leg3->tuple.depth = leg2.tuple.depth + positionshift;
      store ( leg3->occurrences, *legoccurrencesptr );
    }
  }

  if ( (legoccurrencesptr = join ( leg.occurrences )) ) {
    PathLegPtr leg3 = new PathLeg;
    legs.push_back ( leg3 );
    leg3->tuple.connectingnode = leg.tuple.connectingnode;
    leg3->tuple.edgelabel = leg.tuple.edgelabel;
    leg3->tuple.nodelabel = leg.tuple.nodelabel;
    leg3->tuple.depth = leg.tuple.depth + positionshift;
    store ( leg3->occurrences, *legoccurrencesptr );
  }

  for ( i++; i < parentpath.legs.size (); i++ ) {
    PathLeg &leg2 = (*parentpath.legs[i]);
    if ( (legoccurrencesptr = join ( leg.occurrences, leg2.tuple.connectingnode, leg2.occurrences )) ) {
      PathLegPtr leg3 = new PathLeg;
      legs.push_back ( leg3 );
      leg3->tuple.connectingnode = leg2.tuple.connectingnode;
      leg3->tuple.edgelabel = leg2.tuple.edgelabel;
      leg3->tuple.nodelabel = leg2.tuple.nodelabel;
      leg3->tuple.depth = leg2.tuple.depth + positionshift;
      store ( leg3->occurrences, *legoccurrencesptr );
    }
  }

  if ( positionshift ) {
    addCloseExtensions ( closelegs, leg.occurrences.number ); // stored separately
    return;
  }

  extend ( leg.occurrences );
  for ( unsigned int i = 0; i < fm::candidatelegsoccurrences.size (); i++ ) {
    if ( fm::candidatelegsoccurrences[i].frequency >= fm::minfreq ) {
      PathLegPtr leg2 = new PathLeg;
      legs.push_back ( leg2 );
      leg2->tuple.edgelabel = i;
      leg2->tuple.connectingnode = fm::graphstate->lastNode ();
      DatabaseEdgeLabel &databaseedgelabel = fm::database->edgelabels[fm::database->edgelabelsindexes[i]];
      if ( databaseedgelabel.fromnodelabel == leg.tuple.nodelabel )
        leg2->tuple.nodelabel = databaseedgelabel.tonodelabel;
      else
        leg2->tuple.nodelabel = databaseedgelabel.fromnodelabel;
      leg2->tuple.depth = leg.tuple.depth + 1;
      store ( leg2->occurrences, fm::candidatelegsoccurrences[i] ); // avoid copying
    }
  }

  addCloseExtensions ( closelegs, leg.occurrences.number );
}

Path::~Path () {
  for ( unsigned int i = 0; i < legs.size (); i++ )
    delete legs[i];
  for ( unsigned int i = 0; i < closelegs.size (); i++ )
    delete closelegs[i];
}

// ADDED
bool Path::is_normal ( EdgeLabel edgelabel ) {
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









GSWalk* Path::expand2 (pair<float,string> max, const int parent_size) {

  assert(parent_size>0);

  fm::statistics->patternsize++;
  if ( (unsigned) fm::statistics->patternsize > fm::statistics->frequenttreenumbers.size () ) {
    fm::statistics->frequenttreenumbers.push_back ( 0 );
    fm::statistics->frequentpathnumbers.push_back ( 0 );
    fm::statistics->frequentgraphnumbers.push_back ( 0 );
  }
  ++fm::statistics->frequentpathnumbers[fm::statistics->patternsize-1];
  
  if ( fm::statistics->patternsize == ((1<<(sizeof(NodeId)*8))-1) ) {
    fm::statistics->patternsize--;
    return new GSWalk();
  }

  vector<unsigned int> forwpathlegs; forwpathlegs.clear();
  vector<unsigned int> backwpathlegs; backwpathlegs.clear();
  vector<unsigned int> pathlegs; pathlegs.clear();

  for ( unsigned int i = 0; i < legs.size (); i++ ) {

    PathTuple &tuple = legs[i]->tuple;
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
    PathTuple &tuple = legs[i]->tuple;
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
 

  // horizontal view: conflict_resolution will merge into siblingwalk
  // NOTE: siblingwalk is intended to 'carry' the growing meta pattern
  GSWalk* siblingwalk = new GSWalk();

  vector<int> core_ids; 
  for (int j=0; j<parent_size; j++) core_ids.push_back(j);
  int legcnt=0;
  
  // Grow Path forw
  for (unsigned int j=0; j<forwpathlegs.size() ; j++ ) {
    unsigned int index = forwpathlegs[j];

    GSWalk* gsw = new GSWalk();
    GSWalk* topdown = NULL;

    bool nsign=1;

    #ifdef DEBUG
    int diehard = 0;
    #endif

    // Calculate chisq
    if (fm::chisq->active) fm::chisq->Calc(legs[index]->occurrences.elements);
    float cur_chisq=fm::chisq->p;
          
    // GRAPHSTATE AND OUTPUT
    fm::graphstate->insertNode ( legs[index]->tuple.connectingnode, legs[index]->tuple.edgelabel, legs[index]->occurrences.maxdegree );
    #ifdef DEBUG
    fm::graphstate->print(legs[index]->occurrences.frequency);
    #endif

    // immediate output

    #ifdef DEBUG
    fm::gsp_out=false;
    string s = fm::graphstate->to_s(legs[index]->occurrences.frequency);
    if (s.find("C-C=C-O-C-N")!=string::npos) { fm::die=1; diehard=1; }
    fm::die=1;
    #endif
   
    if (fm::chisq->active) {
        map<Tid, int> weightmap_a; each_it(fm::chisq->fa_set, set<Tid>::iterator) { weightmap_a.insert(make_pair((*it),1)); }
        map<Tid, int> weightmap_i; each_it(fm::chisq->fi_set, set<Tid>::iterator) { weightmap_i.insert(make_pair((*it),1)); }
        fm::graphstate->print(gsw, weightmap_a, weightmap_i);
        gsw->activating=fm::chisq->activating;
        if (cur_chisq >= fm::chisq->sig) {
            nsign=0;
        }
    }
    const int gsw_size=gsw->nodewalk.size();

    // !STOP: MERGE TO SIBLINGWALK
    if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Already nodes marked as available 2.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }

    if (nsign || gsw->activating!=siblingwalk->activating) {
          if (siblingwalk->hops>0) {
              if (siblingwalk->hops>1) {
                  siblingwalk->svd();
              }
              cout << siblingwalk;
          }
          delete siblingwalk;
          siblingwalk = new GSWalk();
    }
    if (!nsign && ((gsw->activating==siblingwalk->activating) || !siblingwalk->edgewalk.size())) {
        #ifdef DEBUG
        if (fm::die) cout << "CR gsw 1" << endl;
        #endif
        int res=gsw->conflict_resolution(core_ids, siblingwalk);
    }

    if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Still nodes marked as available 2.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }

    // RECURSE
    if ( ( !fm::do_pruning || (fm::chisq->u >= fm::chisq->sig) ) &&
         (  fm::refine_singles || (legs[index]->occurrences.frequency>1) )
       ) {   // UB-PRUNING
            Path path ( *this, index );
            if (max.first<fm::chisq->p) { fm::updated = true; topdown = path.expand2 ( pair<float, string>(fm::chisq->p, fm::graphstate->to_s(legs[index]->occurrences.frequency)), gsw_size); }
            else topdown = path.expand2 (max,  gsw_size);
    }

    // merge to siblingwalk
    if (topdown != NULL) {
         if (topdown->edgewalk.size()) {

              #ifdef DEBUG
              if (fm::die) {
                  cout << "TOPDOWN BEGIN " << core_ids.size() << endl;
                  cout << topdown ;
                  cout << "--result--" << endl;
                  cout << siblingwalk ;
              }
              #endif

              if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Already nodes marked as available 2.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }
              // STOP: OUTPUT TOPDOWN
              if (nsign || siblingwalk->activating!=topdown->activating) {
                  #ifdef DEBUG
                  if (fm::die) cout << "STOP CRITERIUM at CHI " << cur_chisq << endl;
                  #endif
                  if (topdown->hops>0) { 
                      if (topdown->hops>1) { 
                          topdown->svd();
                      }
                      cout << topdown; 
                  }
              }
              // ELSE: MERGE TO SIBLINGWALK
              else {
                  int res=topdown->conflict_resolution(core_ids, siblingwalk);
              }
              if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Still nodes marked as available 2.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; each_it(topdown->to_nodes_ex, vector<int>::iterator) cout << *it << " "; cout << endl;  exit(1); }

              #ifdef DEBUG
              if (fm::die) {
                  cout << "TOPDOWN END " << core_ids.size() << endl;
                  cout << topdown ;
                  cout << "--result--" << endl;
                  cout << siblingwalk ;
              }
              #endif
         }
    }

    fm::graphstate->deleteNode ();

    delete topdown;
    delete gsw;

    #ifdef DEBUG
    if (diehard==1) { 
       cerr << "DYING HARD!" << endl;
       exit(0);
    }
    #endif


  }



  // Grow Path backw
  for (unsigned int j=0; j<backwpathlegs.size() ; j++ ) {
    unsigned int index = backwpathlegs[j];
    
    GSWalk* gsw = new GSWalk();
    GSWalk* topdown = NULL;

    bool nsign=1;

    // Calculate chisq
    if (fm::chisq->active) fm::chisq->Calc(legs[index]->occurrences.elements);
    float cur_chisq = fm::chisq->p;

    // GRAPHSTATE AND OUTPUT
    fm::graphstate->insertNode ( legs[index]->tuple.connectingnode, legs[index]->tuple.edgelabel, legs[index]->occurrences.maxdegree );
    #ifdef DEBUG
    fm::graphstate->print(legs[index]->occurrences.frequency);
    #endif

    if (fm::chisq->active) {
        map<Tid, int> weightmap_a; each_it(fm::chisq->fa_set, set<Tid>::iterator) { weightmap_a.insert(make_pair((*it),1)); }
        map<Tid, int> weightmap_i; each_it(fm::chisq->fi_set, set<Tid>::iterator) { weightmap_i.insert(make_pair((*it),1)); }
        fm::graphstate->print(gsw, weightmap_a, weightmap_i);
        gsw->activating=fm::chisq->activating;
        if (cur_chisq >= fm::chisq->sig) {
            nsign=0;
        }
    }
    const int gsw_size=gsw->nodewalk.size();

    // !STOP: MERGE TO SIBLINGWALK
    if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Already nodes marked as available 3.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }

    if (nsign || gsw->activating!=siblingwalk->activating) {
          if (siblingwalk->hops>0) {
              if (siblingwalk->hops>1) {
                  siblingwalk->svd();
              }
              cout << siblingwalk;
          }
          delete siblingwalk;
          siblingwalk = new GSWalk();
    }
    if (!nsign && ((gsw->activating==siblingwalk->activating) || !siblingwalk->edgewalk.size())) {
        #ifdef DEBUG
        if (fm::die) cout << "CR gsw 2" << endl;
        #endif
        int res=gsw->conflict_resolution(core_ids, siblingwalk);
    }

    if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Still nodes marked as available 3.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }
 

    // RECURSE
    if ( ( !fm::do_pruning || (fm::chisq->u >= fm::chisq->sig) ) &&
         (  fm::refine_singles || (legs[index]->occurrences.frequency>1) )
       ) {   // UB-PRUNING
            Path path ( *this, index );
            if (max.first<fm::chisq->p) { fm::updated = true; topdown = path.expand2 ( pair<float, string>(fm::chisq->p, fm::graphstate->to_s(legs[index]->occurrences.frequency)), gsw_size); }
            else topdown = path.expand2 (max, gsw_size);
    }

    // merge to siblingwalk
    if (topdown != NULL) {
         if (topdown->edgewalk.size()) {

              #ifdef DEBUG
              if (fm::die) {
                  cout << "TOPDOWN BEGIN " << core_ids.size() << endl;
                  cout << topdown ;
                  cout << "--result--" << endl;
                  cout << siblingwalk ;
              }
              #endif

              if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Already nodes marked as available 3.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }
              // STOP: OUTPUT TOPDOWN
              if (nsign || siblingwalk->activating!=topdown->activating) {
                  #ifdef DEBUG
                  if (fm::die) cout << "STOP CRITERIUM at CHI " << cur_chisq << endl;
                  #endif
                  if (topdown->hops>0) { 
                      if (topdown->hops>1) { 
                          topdown->svd();
                      }
                      cout << topdown; 
                  }
              }
              // ELSE: MERGE TO SIBLINGWALK
              else {
                  int res=topdown->conflict_resolution(core_ids, siblingwalk); 
              }
              if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Still nodes marked as available 3.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }

              #ifdef DEBUG
              if (fm::die) {
                  cout << "TOPDOWN END " << core_ids.size() << endl;
                  cout << topdown ;
                  cout << "--result--" << endl;
                  cout << siblingwalk ;
              }
              #endif
         }
    }

   
    fm::graphstate->deleteNode ();

    delete topdown;
    delete gsw;

  }


  bool uptmp = fm::updated;

  if (fm::bbrc_sep && legs.size() > 0) {
      if (fm::do_output && !fm::console_out && fm::result->size() && (fm::result->back()!=fm::graphstate->sep())) (*fm::result) << fm::graphstate->sep();
  }

  for ( unsigned int i = 0; i < legs.size (); i++ ) {
    PathTuple &tuple = legs[i]->tuple;
    if ( tuple.depth != nodelabels.size () - 1 ) {

    // PHASE 2: GROW TREE
    if ( legs[i]->tuple.depth != 0 ) {
      if ( ( totalsymmetry || legs[i]->tuple.depth <= edgelabels.size () / 2 ) &&
      ( legs[i]->tuple.depth != 1 || legs[i]->tuple.edgelabel >= edgelabels[0] ) &&
      ( legs[i]->tuple.depth != nodelabels.size () - 2 || legs[i]->tuple.edgelabel >= edgelabels.back () ) &&
      fm::type > 1 ) {

          // new current pattern
          GSWalk* gsw = new GSWalk();
          GSWalk* topdown = NULL;

          bool nsign=1;

          if (fm::chisq->active) fm::chisq->Calc(legs[i]->occurrences.elements);
          float cur_chisq = fm::chisq->p;

          fm::graphstate->insertNode ( legs[i]->tuple.connectingnode, legs[i]->tuple.edgelabel, legs[i]->occurrences.maxdegree );
          #ifdef DEBUG
          fm::graphstate->print(legs[i]->occurrences.frequency);
          #endif

          #ifdef DEBUG
          fm::gsp_out=false;
          string s = fm::graphstate->to_s(legs[i]->occurrences.frequency);
          bool diehard=0;
          //if (s.find("C-C(-O-C-N-O)(=C-C)")!=string::npos) { fm::die=1; diehard=1; }
          #endif

          if (fm::chisq->active) {
              map<Tid, int> weightmap_a; each_it(fm::chisq->fa_set, set<Tid>::iterator) { weightmap_a.insert(make_pair((*it),1)); }
              map<Tid, int> weightmap_i; each_it(fm::chisq->fi_set, set<Tid>::iterator) { weightmap_i.insert(make_pair((*it),1)); }
              fm::graphstate->print(gsw, weightmap_a, weightmap_i);
              gsw->activating=fm::chisq->activating;
              if (cur_chisq >= fm::chisq->sig) {
                  nsign=0;
              }
          }
          const int gsw_size=gsw->nodewalk.size();

          // !STOP: MERGE TO SIBLINGWALK
          if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Already nodes marked as available 4.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }

          if (nsign || gsw->activating!=siblingwalk->activating) {
                if (siblingwalk->hops>0) {
                    if (siblingwalk->hops>1) {
                        siblingwalk->svd();
                    }
                    cout << siblingwalk;
                }
                delete siblingwalk;
                siblingwalk = new GSWalk();
          }
          if (!nsign && ((gsw->activating==siblingwalk->activating) || !siblingwalk->edgewalk.size())) {
              #ifdef DEBUG
              if (fm::die) cout << "CR gsw 3" << endl;
              #endif
              int res=gsw->conflict_resolution(core_ids, siblingwalk);
          }


          if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Still nodes marked as available 4.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }

          if ( ( !fm::do_pruning ||  (fm::chisq->u >= fm::chisq->sig ) ) &&
               (  fm::refine_singles || (legs[i]->occurrences.frequency>1) )
             ) {
              PatternTree tree ( *this, i );
              if (max.first<cur_chisq) { fm::updated = true; topdown = tree.expand ( pair<float, string>(cur_chisq, fm::graphstate->to_s(legs[i]->occurrences.frequency)), gsw_size); }
              else topdown = tree.expand (max, gsw_size);
          }

          // merge to siblingwalk
          if (topdown != NULL) {
               if (topdown->edgewalk.size()) {

                    #ifdef DEBUG
                    if (fm::die) {
                        cout << "TOPDOWN BEGIN " << core_ids.size() << endl;
                        cout << topdown ;
                        cout << "--result--" << endl;
                        cout << siblingwalk ;
                    }
                    #endif

                    if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Already nodes marked as available 4.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }
                    // STOP: OUTPUT TOPDOWN
                    if (nsign || siblingwalk->activating!=topdown->activating) {
                        #ifdef DEBUG
                        if (fm::die) cout << "STOP CRITERIUM at CHI " << cur_chisq << endl;
                        #endif
                        if (topdown->hops>0) { 
                            if (topdown->hops>1) { 
                                topdown->svd();
                            }
                            cout << topdown; 
                        }
                    }
                    // ELSE: MERGE TO SIBLINGWALK
                    else {
                        int res=topdown->conflict_resolution(core_ids, siblingwalk);
                    }
                    if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Still nodes marked as available 4.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }

                    #ifdef DEBUG
                    if (fm::die) {
                        cout << "TOPDOWN END " << core_ids.size() << endl;
                        cout << topdown;
                        cout << "--result--" << endl;
                        cout << siblingwalk ;
                    }
                    #endif
               }
          }


	      fm::graphstate->deleteNode ();
          delete topdown;
          delete gsw;
          #ifdef DEBUG
          if (diehard==1) { 
             cerr << "DYING HARD! " << legs.size() << endl;
             exit(0);
          }
          #endif

        }

        #ifdef DEBUG  
        if (!legs.size()) cout << fm::graphstate->sep() << endl;
        #endif
      }
    }
  }

  // delete horizontal view
  fm::updated=uptmp;
  fm::statistics->patternsize--;
  return siblingwalk;

//  cerr << "backtracking p" << endl;
}








void Path::expand () {

  //fm::die=1;
  // horizontal view: conflict_resolution will merge into siblingwalk
  // NOTE: siblingwalk is intended to 'carry' the growing meta pattern
  GSWalk* siblingwalk = new GSWalk();
  vector<int> core_ids; core_ids.push_back(0); core_ids.push_back(1);
  int legcnt=0;

  for ( unsigned int i = 0; i < legs.size (); i++ ) {

    GSWalk* gsw = new GSWalk(); 
    GSWalk* topdown = NULL;

    bool nsign=1;

    PathTuple &tuple = legs[i]->tuple;
    if ( tuple.nodelabel >= nodelabels[0] ) {
        
      if (fm::chisq->active) fm::chisq->Calc(legs[i]->occurrences.elements);
      float cur_chisq = fm::chisq->p;

      // GRAPHSTATE AND OUTPUT
      fm::graphstate->insertNode ( tuple.connectingnode, tuple.edgelabel, legs[i]->occurrences.maxdegree );
      #ifdef DEBUG
      fm::graphstate->print(legs[i]->occurrences.frequency);
      #endif

      if (fm::chisq->active) {
          map<Tid, int> weightmap_a; each_it(fm::chisq->fa_set, set<Tid>::iterator) { weightmap_a.insert(make_pair((*it),1)); }
          map<Tid, int> weightmap_i; each_it(fm::chisq->fi_set, set<Tid>::iterator) { weightmap_i.insert(make_pair((*it),1)); }
          fm::graphstate->print(gsw, weightmap_a, weightmap_i);
          gsw->activating=fm::chisq->activating;
          if (cur_chisq >= fm::chisq->sig) {
              nsign=0;
          }
      }
      const int gsw_size=gsw->nodewalk.size();

      // !STOP: MERGE TO SIBLINGWALK
      if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Already nodes marked as available 1.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }

      if (nsign || gsw->activating!=siblingwalk->activating) {
            if (siblingwalk->hops>0) {
                if (siblingwalk->hops>1) {
                    siblingwalk->svd();
                }
                cout << siblingwalk;
            }
            delete siblingwalk;
            siblingwalk = new GSWalk();
      }
      if (!nsign && ((gsw->activating==siblingwalk->activating) || !siblingwalk->edgewalk.size())) {
          #ifdef DEBUG
          if (fm::die) cout << "CR gsw 4" << endl;
          #endif
          int res=gsw->conflict_resolution(core_ids, siblingwalk);
      }

      if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Still nodes marked as available 1.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }


      // RECURSE
      Path path (*this, i);
      fm::updated = true;

      topdown = path.expand2 (pair<float, string>(fm::chisq->p, fm::graphstate->to_s(legs[i]->occurrences.frequency)), gsw_size);

      // merge to siblingwalk
      if (topdown != NULL) {
           if (topdown->edgewalk.size()) {

                #ifdef DEBUG
                if (fm::die) {
                    cout << "TOPDOWN BEGIN " << core_ids.size() << endl;
                    cout << topdown;
                    cout << "--result--" << endl;
                    cout << siblingwalk ;
                }
                #endif

                if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Already nodes marked as available 1.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }
                // STOP: OUTPUT TOPDOWN
                if (nsign || siblingwalk->activating!=topdown->activating) {
                    #ifdef DEBUG
                    if (fm::die) cout << "STOP CRITERIUM at CHI " << cur_chisq << endl;
                    #endif
                    if (topdown->hops>0) { 
                        if (topdown->hops>1) { 
                            topdown->svd();
                        }
                        cout << topdown; 
                    }
                }
                // ELSE: MERGE TO SIBLINGWALK
                else {
                    int res=topdown->conflict_resolution(core_ids, siblingwalk); 
                }
                if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Still nodes marked as available 1.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }

                #ifdef DEBUG
                if (fm::die) {
                    cout << "TOPDOWN END " << core_ids.size() << endl;
                    cout << topdown ;
                    cout << "--result--" << endl;
                    cout << siblingwalk ;
                }
                #endif
           }
      }

      fm::graphstate->deleteNode ();

    }

    delete gsw;    
    delete topdown;

  }
  fm::graphstate->deleteStartNode ();
  delete siblingwalk;

//  cerr << "backtracking p" << endl;
}






ostream &operator<< ( ostream &stream, Path &path ) {
  stream << /* database->nodelabels[ */ (int) path.nodelabels[0] /* ].inputlabel; */ << " ";
  for ( unsigned int i = 0; i < path.edgelabels.size (); i++ ) {
    //stream << (char) ( path.edgelabels[i] + 'A' ) << path.nodelabels[i+1];
    stream << /*database->edgelabels[database->edgelabelsindexes[*/ (int) path.edgelabels[i] /*]].inputedgelabel */ << " " <<  /* database->nodelabels[ */ (int) path.nodelabels[i+1] /* ].inputlabel */ << " ";
  }
  return stream;
}
