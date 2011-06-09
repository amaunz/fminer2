// patterntree.cpp
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

#include "patterntree.h"
#include "graphstate.h"

namespace fm {
    extern unsigned int last_minfreq;
    extern bool last_updated;
    extern bool last_do_pruning;
    extern bool last_console_out;
    extern bool last_refine_singles;
    extern bool last_do_output;
    extern bool last_bbrc_sep;
    extern bool last_regression;
    extern bool last_gsp_out;
    extern int last_die;
    extern bool last_do_last;

    extern LastDatabase* last_database;
    extern ChisqLastConstraint* last_chisq;
    extern KSLastConstraint* last_ks;
    extern vector<string>* last_result;
    extern LastStatistics* last_statistics;
    extern LastGraphState* last_graphstate;
    extern LastLegOccurrences* last_legoccurrences;

    extern vector<LastLegOccurrences> last_Lastcandidatelegsoccurrences; 
    extern int last_max_hops;
}

int maxsize = ( 1 << ( sizeof(LastNodeId)*8 ) ) - 1; // safe default for the largest allowed pattern

inline void LastPatternTree::addLastLeg ( LastNodeId connectingnode, const int depth, const LastEdgeLabel edgelabel, LastLegOccurrences &legoccurrences ) {
  LastLegPtr leg = new LastLeg;
  leg->tuple.depth = depth;
  leg->tuple.label = edgelabel;
  leg->tuple.connectingnode = connectingnode;
  store ( leg->occurrences, legoccurrences );
  legs.push_back ( leg );
}

// this function assumes that the extension tuple is already added at the back of the queue,
// and the equivalency information has been filled in.
void LastPatternTree::addExtensionLastLegs ( LastTuple &tuple, LastLegOccurrences &legoccurrences ) {
  if ( legoccurrences.maxdegree == 1 )
    return;
  if ( tuple.depth == maxdepth ) {
    bbrc_extend ( legoccurrences, MAXEDGELABEL, (unsigned char) NONODE );
    LastaddCloseExtensions ( closelegs, legoccurrences.number );
    return;
  }
  LastEdgeLabel minlabel = NOEDGELABEL, neglect = '\0', pathlowestlabel = treetuples[tuple.depth + 1 + rootpathstart].label;

  if ( nextprefixindex != NONEXTPREFIX ) {
    if ( treetuples[nextprefixindex].depth <= tuple.depth ) {
      // heuristic saving
      bbrc_extend ( legoccurrences, MAXEDGELABEL, (unsigned char) NONODE );
      LastaddCloseExtensions ( closelegs, legoccurrences.number );
      return;
    }
    minlabel = treetuples[nextprefixindex].label;
    if ( minlabel == pathlowestlabel )
      minlabel = NOEDGELABEL;
    else
      neglect = pathlowestlabel;
  }
  if ( tuple.depth == maxdepth - 1 ) {
    if ( rootpathrelations.back () > 0 ) {
      // heuristic saving
      bbrc_extend ( legoccurrences, MAXEDGELABEL, (unsigned char) NONODE );
      LastaddCloseExtensions ( closelegs, legoccurrences.number );
      return;
    }
    if ( rootpathrelations.back () == 0 )
      if ( minlabel != NOEDGELABEL ) {
        neglect = NOEDGELABEL;
    	if ( pathlowestlabel > minlabel )
          minlabel = pathlowestlabel;
      }
      else {
        neglect = NOEDGELABEL;
        minlabel = pathlowestlabel;
      }
      // else we have a restriction, and that restriction is at least as high
      // as the label on the path.
  }

  if ( minlabel != NOEDGELABEL )
    bbrc_extend ( legoccurrences, minlabel, neglect );
  else
    bbrc_extend ( legoccurrences );

  if ( fm::last_Lastcandidatelegsoccurrences[pathlowestlabel].frequency >= fm::last_minfreq )
    // this is the first possible extension, as we force this label to be the lowest!
    addLastLeg ( fm::last_graphstate->lastNode (), tuple.depth + 1, pathlowestlabel, fm::last_Lastcandidatelegsoccurrences[pathlowestlabel] );

  for ( int i = 0; (unsigned) i < fm::last_Lastcandidatelegsoccurrences.size (); i++ ) {
    if ( fm::last_Lastcandidatelegsoccurrences[i].frequency >= fm::last_minfreq && i != pathlowestlabel )
      addLastLeg ( fm::last_graphstate->lastNode (), tuple.depth + 1, i, fm::last_Lastcandidatelegsoccurrences[i] );
  }

  LastaddCloseExtensions ( closelegs, legoccurrences.number );
}

void LastPatternTree::addLeftLastLegs ( LastPath &path, LastPathLastLeg &leg, int &i, LastDepth olddepth, LastEdgeLabel lowestlabel, int leftend, int edgesize2 ) {
  // the order of the other extensions is almost correct - except that we must move
  // extensions at depth d which are for the label of the path at depth d to the front of
  // all legs of depth d...

  for ( ; i < (int) path.legs.size () && (int) path.legs[i]->tuple.depth <= leftend; i++ ) {
    if ( path.legs[i]->tuple.depth != olddepth ) {
      // when we encounter a new depth, we look whether there is an extension tuple at that
      // depth with the same label as the node on the path; this tuple has to be moved to the front; ...
      olddepth = path.legs[i]->tuple.depth;
      lowestlabel = path.edgelabels[olddepth - 1];
      int i2 = i;
      while ( (unsigned) i2 < path.legs.size () && path.legs[i2]->tuple.depth == olddepth ) {
        if ( path.legs[i2]->tuple.edgelabel == lowestlabel ) {
          LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[i2]->tuple.connectingnode, path.legs[i2]->occurrences );
          if ( legoccurrencesptr )
            addLastLeg ( path.legs[i2]->tuple.connectingnode, edgesize2 - olddepth, path.legs[i2]->tuple.edgelabel, *legoccurrencesptr );
          break;
        }
        i2++;
      }
    }
    // skip lowest label tuples, as they have already been moved to the front...
    if ( path.legs[i]->tuple.edgelabel != lowestlabel ) {
      LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[i]->tuple.connectingnode, path.legs[i]->occurrences );
      if ( legoccurrencesptr )
        addLastLeg ( path.legs[i]->tuple.connectingnode, edgesize2 - olddepth, path.legs[i]->tuple.edgelabel, *legoccurrencesptr );
    }
  }
}

int LastPatternTree::addLeftLastLegs ( LastPath &path, LastPathLastLeg &leg, LastTuple &tuple, unsigned int legindex, int leftend, int edgesize2 ) {
  int i;

  LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences );
  if ( legoccurrencesptr )
    addLastLeg ( leg.tuple.connectingnode, tuple.depth, tuple.label, *legoccurrencesptr );

  // the easy part - the extensions of the left side of the original path

  if ( rootpathrelations.back () == 0 && tuple.depth != maxdepth ) {
    // we force the label on the path to be the lowest label, so all labels which are actually
    // lower must be added as they are higher here, deepest nodes are an exception;
    // no path may be lower than the current path

    for ( i = legindex - 1; i >= 0 && path.legs[i]->tuple.depth == leg.tuple.depth; i-- );
    LastEdgeLabel lowestlabel = path.edgelabels[leg.tuple.depth - 1];
    for ( i++; i < (int) legindex; i++ ) {
      if ( path.legs[i]->tuple.edgelabel != lowestlabel ) {
        legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[i]->tuple.connectingnode, path.legs[i]->occurrences );
        if ( legoccurrencesptr )
          addLastLeg ( path.legs[i]->tuple.connectingnode, tuple.depth, path.legs[i]->tuple.edgelabel, *legoccurrencesptr );
      }
    }
  }

  // the order of the other extensions is almost correct - except that we must move
  // extensions at depth d which are for the label of the path at depth d to the front of
  // all legs of depth d...

  i = legindex + 1;
  addLeftLastLegs ( path, leg, i, leg.tuple.depth, path.edgelabels[leg.tuple.depth - 1], leftend, edgesize2 );

  return i;
}

void LastPatternTree::addRightLastLegs ( LastPath &path, LastPathLastLeg &leg, int &i, LastDepth olddepth, LastEdgeLabel lowestlabel, int rightstart, int nodesize2 ) {
  int i2 = i + 1, k;
  LastLegOccurrencesPtr legoccurrencesptr;

  while ( i >= 0 && (int) path.legs[i]->tuple.depth >= rightstart ) {
    // we encounter a new depth, do all the extensions of the previous depths...
    if ( path.legs[i]->tuple.depth != olddepth ) {
      for ( k = i + 1; k < i2; k++ ) {
        if ( path.legs[k]->tuple.edgelabel != lowestlabel ) {
          legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
          if ( legoccurrencesptr )
            addLastLeg ( path.legs[k]->tuple.connectingnode, path.legs[k]->tuple.depth - nodesize2, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
        }
      }
      i2 = i + 1;
      olddepth = path.legs[i]->tuple.depth;
      lowestlabel = path.edgelabels[olddepth];
    }
    // for the current depth, we encounter the label that is on the second path,
    // this extension must be moved to the front, so before the other tuples are
    // added by the code above
    if ( path.legs[i]->tuple.edgelabel == lowestlabel ) {
      legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[i]->tuple.connectingnode, path.legs[i]->occurrences );
      if ( legoccurrencesptr )
        addLastLeg ( path.legs[i]->tuple.connectingnode, path.legs[i]->tuple.depth - nodesize2, path.legs[i]->tuple.edgelabel, *legoccurrencesptr );
    }
    i--;
  }
  // some tuples may not have been checked yet
  for ( k = i + 1; k < i2; k++ ) {
    if ( path.legs[k]->tuple.edgelabel != lowestlabel ) {
      legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
      if ( legoccurrencesptr )
        addLastLeg ( path.legs[k]->tuple.connectingnode, path.legs[k]->tuple.depth - nodesize2, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
    }
  }
}

int LastPatternTree::addRightLastLegs ( LastPath &path, LastPathLastLeg &leg, LastTuple &tuple, unsigned int legindex, int rightstart, int nodesize2 ) {
  int i;
  LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences );
  if ( legoccurrencesptr )
    addLastLeg ( leg.tuple.connectingnode, tuple.depth, tuple.label, *legoccurrencesptr );

  // other extensions at the right path

  i = legindex - 1;
  while ( i >= 0 && path.legs[i]->tuple.depth == leg.tuple.depth )
    i--;

  // brothers of the new leg
  if ( rootpathrelations.back () == 0 && tuple.depth != maxdepth )
    for ( int j = i + 1; j < (int) legindex; j++ ) {
      LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[j]->tuple.connectingnode, path.legs[j]->occurrences );
      if ( legoccurrencesptr )
        addLastLeg ( path.legs[j]->tuple.connectingnode, tuple.depth, path.legs[j]->tuple.edgelabel, *legoccurrencesptr );
    }
  LastEdgeLabel lowestlabel = path.edgelabels[leg.tuple.depth];
  for ( int j = legindex + 1; j < (int) path.legs.size () && path.legs[j]->tuple.depth == leg.tuple.depth; j++ ) {
    if ( path.legs[j]->tuple.edgelabel != lowestlabel ) {
      LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[j]->tuple.connectingnode, path.legs[j]->occurrences );
      if ( legoccurrencesptr )
        addLastLeg ( path.legs[j]->tuple.connectingnode, tuple.depth, path.legs[j]->tuple.edgelabel, *legoccurrencesptr );
    }
  }

  // candidates at a lower depth
  if ( i >= 0 )
    addRightLastLegs ( path, leg, i, leg.tuple.depth, lowestlabel, rightstart, nodesize2 );

  return i;
}

LastPatternTree::LastPatternTree ( LastPath &path, unsigned int legindex ) {
  LastPathLastLeg &leg = (*path.legs[legindex]);
  
  maxdepth = path.edgelabels.size () / 2 - 1;
  int leftwalk, leftstart, rightwalk, rightstart;
  LastLegOccurrencesPtr legoccurrencesptr;

  LastaddCloseExtensions ( closelegs, path.closelegs, leg.occurrences );

  int nodesize2 = path.nodelabels.size () / 2;
  int edgesize2 = path.edgelabels.size () / 2;
  rightstart = rightwalk = nodesize2;
  leftstart = leftwalk = edgesize2 - 1;
  if ( path.totalsymmetry == 0 )
    symmetric = path.nodelabels.size () % 2 + 1;
  else
    symmetric = 0;
  if ( path.totalsymmetry ||
       leg.tuple.depth * 2 == path.edgelabels.size () ) {
    while ( leftwalk >= 0 &&
            path.edgelabels[leftwalk] == path.edgelabels[rightwalk] )
      leftwalk--, rightwalk++;
      // now there is one very nasty case: A-B-A-B-A-B
      // this path is not symmetric, but the labels on the edges are!
      // So we can find that leftwalk == -1.
      // In this case, we assume that the left part is the first path,
      // furthermore the position of the extension determines to which path
      // it is added
    fm::last_graphstate->nasty = ( leftwalk == -1 );
    if ( leftwalk == -1 || path.edgelabels[leftwalk] < path.edgelabels[rightwalk] ) {
      // left part of the path should be the first path in the tree

      treetuples.reserve ( path.edgelabels.size () + 1 );
      int i, j;
      for ( i = leftstart, j = 0; i >= 0; i--, j++ ) {
        vector_push_back ( LastTuple, treetuples, tuple );
        tuple.label = path.edgelabels[i];
        tuple.depth = j;
      }

      int leftend = nodesize2 - 1;

      if ( (int) leg.tuple.depth <= leftend ) {
        // extension is an extension of the first path

        // add extension tuple
        vector_push_back ( LastTuple, treetuples, tuple );
        tuple.label = leg.tuple.edgelabel;
        tuple.depth = edgesize2 - leg.tuple.depth;

        rootpathrelations.reserve ( tuple.depth );
        rightmostindexes.reserve ( tuple.depth );
        for ( i = 0; i < (int) tuple.depth; i++ ) {
          rootpathrelations.push_back ( 0 );
          rightmostindexes.push_back ( i );
        }
        rightmostindexes.push_back ( treetuples.size () - 1 );
        if ( tuple.label == treetuples[tuple.depth].label ) {
          nextprefixindex = tuple.depth + 1;
          rootpathrelations.push_back ( 0 );
        }
        else {
          nextprefixindex = NONEXTPREFIX;
          rootpathrelations.push_back ( treetuples[tuple.depth].label - tuple.label );
        }

        rootpathstart = 0;
        nextpathstart = treetuples.size ();

        // fill in possible extensions below the new node
        addExtensionLastLegs ( tuple, leg.occurrences );

        // add second path
        for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
          vector_push_back ( LastTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }

        // here, we have determined the expansion tuple sequence of the tree, and the normal
        // form information. Now we have to fill in the legs at other positions than below
        // the new node. This is complicated as the occurrence list order of the path
        // has to be changed...

        i = addLeftLastLegs ( path, leg, tuple, legindex, leftend, edgesize2 );

	secondpathleg = legs.size ();


        // a difficult part - the extensions of the right side of the original path; here the labels
        // are in correct order, but the depth order is incorrect (more precisely, reversed)

        int j2 = path.legs.size ();
        int j = j2 - 1, k;

        if ( j >= i && path.legs[j]->tuple.depth == path.edgelabels.size () ) {
          // do not include legs after the longest path!
          do
            j--;
          while ( j >= i && path.legs[j]->tuple.depth == path.edgelabels.size () );
          j2 = j + 1;
        }


        if ( j >= i && path.legs[j]->tuple.depth == path.nodelabels.size () - 2 ) {
          // at the lowest type of the second path, no edge with a lower label than the last on
          // that path may be added. Incorporate that into computation!
          if ( path.legs[j]->tuple.edgelabel >= path.edgelabels.back () ) {
            do {
              j--;
            }
            while ( j >= i && path.legs[j]->tuple.depth == path.nodelabels.size () - 2
                            && path.legs[j]->tuple.edgelabel >= path.edgelabels.back () );
            for ( k = j + 1; k < j2; k++ ) {
              legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
              if ( legoccurrencesptr ) {
                addLastLeg ( path.legs[k]->tuple.connectingnode, path.legs[k]->tuple.depth - nodesize2, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
              }
            }
          }
          while ( j >= i && path.legs[j]->tuple.depth == path.nodelabels.size () - 2 )
            j--;
        }

        if ( j >= i )
          addRightLastLegs ( path, leg, j, NODEPTH, NOEDGELABEL, rightstart, nodesize2 );

        // we're done!

      }
      else {
        // extension is an extension of the second path

        rootpathstart = nextpathstart = treetuples.size ();

        // add second path
        for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
          vector_push_back ( LastTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }
        // add extension tuple
        vector_push_back ( LastTuple, treetuples, tuple );
        tuple.label = leg.tuple.edgelabel;
        tuple.depth = leg.tuple.depth - nodesize2;

        rootpathrelations.reserve ( tuple.depth + 1 );
        rightmostindexes.reserve ( tuple.depth + 1 );
        for ( i = 0; i < (int) tuple.depth; i++ ) {
          rootpathrelations.push_back ( 0 );
          rightmostindexes.push_back ( i + rootpathstart );
        }
        rightmostindexes.push_back ( treetuples.size () - 1 );
        if ( tuple.label == treetuples[tuple.depth + rootpathstart ].label ) {
          nextprefixindex = tuple.depth + rootpathstart + 1;
          rootpathrelations.push_back ( 0 );
        }
        else {
          nextprefixindex = NONEXTPREFIX;
          rootpathrelations.push_back ( treetuples[tuple.depth + rootpathstart].label - tuple.label );
        }

        // fill in possible extensions below the new node
        addExtensionLastLegs ( tuple, leg.occurrences );

        // add the leg itself

        addRightLastLegs ( path, leg, tuple, legindex, rightstart, nodesize2 );

        secondpathleg = legs.size ();
        // we're done
      }
    }
    else {
      // right part of the tree is the first path in the tree

      treetuples.reserve ( path.edgelabels.size () + 1 );
      int i, j;
      for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
        vector_push_back ( LastTuple, treetuples, tuple );
        tuple.label = path.edgelabels[i];
        tuple.depth = j;
      }

      int leftend = edgesize2;
      rightstart = leftend + 1;

      if ((int) leg.tuple.depth <= leftend ) {
        // extension is an extension of the left, second path
        rootpathstart = nextpathstart = treetuples.size ();

        // add second path (the left path)
        for ( i = leftstart, j = 0; i >= 0; i--, j++ ) {
          vector_push_back ( LastTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }
        // add extension tuple
        vector_push_back ( LastTuple, treetuples, tuple );
        tuple.label = leg.tuple.edgelabel;
        tuple.depth = edgesize2 - leg.tuple.depth;

        rootpathrelations.reserve ( tuple.depth + 1 );
        rightmostindexes.reserve ( tuple.depth + 1 );
        for ( i = 0; i < (int) tuple.depth; i++ ) {
          rootpathrelations.push_back ( 0 );
          rightmostindexes.push_back ( i + rootpathstart );
        }
        rightmostindexes.push_back ( treetuples.size () - 1 );
        if ( tuple.label == treetuples[tuple.depth + rootpathstart ].label ) {
          nextprefixindex = tuple.depth + rootpathstart + 1;
          rootpathrelations.push_back ( 0 );
        }
        else {
          nextprefixindex = NONEXTPREFIX;
          rootpathrelations.push_back ( treetuples[tuple.depth + rootpathstart].label - tuple.label );
        }


        // fill in possible extensions below the new node
        addExtensionLastLegs ( tuple, leg.occurrences );

        // add the leg itself, and all legs on the left path above in the tree

        addLeftLastLegs ( path, leg, tuple, legindex, leftend, edgesize2 );

        secondpathleg = legs.size ();

        // we're done
      }
      else {
        // the right path is the lowest, and we're adding a leg in that path
        // add extension tuple

        int i, j;

        vector_push_back ( LastTuple, treetuples, tuple );
        tuple.label = leg.tuple.edgelabel;
        tuple.depth = leg.tuple.depth - nodesize2;

        rootpathrelations.reserve ( tuple.depth + 1 );
        rightmostindexes.reserve ( tuple.depth + 1 );

        for ( i = 0; i < (int) tuple.depth; i++ ) {
          rootpathrelations.push_back ( 0 );
          rightmostindexes.push_back ( i );
        }
        rightmostindexes.push_back ( treetuples.size () - 1 );
        if ( tuple.label == treetuples[tuple.depth].label ) {
          nextprefixindex = tuple.depth + 1;
          rootpathrelations.push_back ( 0 );
        }
        else {
          nextprefixindex = NONEXTPREFIX;
          rootpathrelations.push_back ( treetuples[tuple.depth].label - tuple.label );
        }

        rootpathstart = 0;
        nextpathstart = treetuples.size ();

        // fill in possible extensions below the new node
        addExtensionLastLegs ( tuple, leg.occurrences );

        // add second path
        for ( i = leftstart, j = 0; i >= 0; i--, j++ ) {
          vector_push_back ( LastTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }

        // here, we have determined the expansion tuple sequence of the tree, and the normal
        // form information. Now we have to fill in the legs at other positions than below
        // the new node. This is complicated as the occurrence list order of the path
        // has to be changed...


        i = addRightLastLegs ( path, leg, tuple, legindex, rightstart, nodesize2 );

        secondpathleg = legs.size ();

        // the easier part - the extensions of the left side of the original path, which is the second
        // path in the tree.

        // no legs at the lowest (highest) type
        j = 0;
        while ( j <= i && path.legs[j]->tuple.depth == 0 )
          j++;

        // siblings of the node at the deepest type must have a higher or equal label to the last
        // label on the path, otherwise the string would become lower
        while ( j <= i && path.legs[j]->tuple.depth == 1
                        && path.legs[j]->tuple.edgelabel < path.edgelabels[0] )
          j++;

        if ( j <= i )
          addLeftLastLegs ( path, leg, j, NODEPTH, NOEDGELABEL, leftend, edgesize2 );

        // we're done!
      }
    }
  }
  else {

    treetuples.reserve ( path.edgelabels.size () + 1 );
    int i, j;
    for ( i = leftstart, j = 0; i >= 0; i--, j++ ) {
      vector_push_back ( LastTuple, treetuples, tuple );
      tuple.label = path.edgelabels[i];
      tuple.depth = j;
    }

    // extension is always an extension of the first path, unless the extension
    // is performed on the middle node of a string of odd number of nodes

    // add extension tuple
    vector_push_back ( LastTuple, treetuples, tuple );
    tuple.label = leg.tuple.edgelabel;
    tuple.depth = edgesize2 - leg.tuple.depth;

    rootpathrelations.reserve ( tuple.depth + 1 );
    rightmostindexes.reserve ( tuple.depth + 1 );
    for ( i = 0; i < (int) tuple.depth; i++ ) {
      rootpathrelations.push_back ( 0 );
      rightmostindexes.push_back ( i );
    }
    rightmostindexes.push_back ( treetuples.size () - 1 );
    if ( tuple.label == treetuples[tuple.depth].label ) {
      nextprefixindex = tuple.depth + 1;
      rootpathrelations.push_back ( 0 );
    }
    else {
      nextprefixindex = NONEXTPREFIX;
      rootpathrelations.push_back ( treetuples[tuple.depth].label - tuple.label );
    }

    rootpathstart = 0;
    nextpathstart = treetuples.size ();

    // fill in possible extensions below the new node
    addExtensionLastLegs ( tuple, leg.occurrences );

    // add second path
    for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
      vector_push_back ( LastTuple, treetuples, tuple );
      tuple.label = path.edgelabels[i];
      tuple.depth = j;
    }

    // here, we have determined the expansion tuple sequence of the tree, and the normal
    // form information. Now we have to fill in the legs at other positions than below
    // the new node. This is complicated as the occurrence list order of the path
    // has to be changed...

    i = addLeftLastLegs ( path, leg, tuple, legindex, rightstart - 1, edgesize2 );

    // a difficult part - the extensions of the right side of the original path; here the labels
    // are in correct order, but the depth order is incorrect (more precisely, reversed)
    // furthermore, as the string is symmetric, the possible legs for the second path are limitted by
    // the first side leg of the first path; the second subtree must be equal, or higher!

    secondpathleg = legs.size ();


    int targetdepth = path.nodelabels.size () - 1 - leg.tuple.depth;

    j = path.legs.size () - 1;

    if ( targetdepth == (int) path.edgelabels.size () - 1 ) {
      while ( j >= i && (int) path.legs[j]->tuple.depth > targetdepth )
        j--;
      int j2 = j;
      while ( j >= i && (int) path.legs[j]->tuple.depth == targetdepth && path.legs[j]->tuple.edgelabel >= leg.tuple.edgelabel )
        j--;
      for ( int k = j + 1; k <= j2; k++ ) {
        LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
        if ( legoccurrencesptr )
          addLastLeg ( path.legs[k]->tuple.connectingnode, tuple.depth, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
      }
      while ( j >= i && (int) path.legs[j]->tuple.depth == targetdepth )
        j--;
    }
    else {
      while ( j >= i && (int) path.legs[j]->tuple.depth > targetdepth )
        j--;
      if ( j >= i ) {
        if ( rootpathrelations.back () != 0 && (int) path.legs[j]->tuple.depth == targetdepth ) {
          int j2 = j;
          LastEdgeLabel lowestlabel = path.edgelabels[targetdepth];
          while ( j >= i &&
	          (int) path.legs[j]->tuple.depth == targetdepth &&
	  	  path.legs[j]->tuple.edgelabel >= tuple.label )
            j--;
          for ( int k = j + 1; k <= j2; k++ )
            if ( path.legs[k]->tuple.edgelabel != lowestlabel ) {
              LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
              if ( legoccurrencesptr )
                addLastLeg ( path.legs[k]->tuple.connectingnode, tuple.depth, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
          }
          while ( j >= i && (int) path.legs[j]->tuple.depth == targetdepth )
            j--;
        }
      }
    }

    if ( j >= i )
      addRightLastLegs ( path, leg, j, NODEPTH, NOEDGELABEL, rightstart, nodesize2 );
        // we're done!

    // symmetry
  }
  
  // ADDED
  fm::last_graphstate->backbonelength = path.nodelabels.size ();
  if ( fm::last_graphstate->backbonelength % 2 == 0 )
    fm::last_graphstate->bicenterlabel = path.edgelabels [ fm::last_graphstate->backbonelength / 2 - 1 ];
  else
    fm::last_graphstate->centerlabel = path.nodelabels [ ( fm::last_graphstate->backbonelength - 1 ) / 2 ];
  fm::last_graphstate->nasty = false;
  fm::last_graphstate->treetuples = &treetuples;
  fm::last_graphstate->closetuples = NULL;
  fm::last_graphstate->startsecondpath = nextpathstart;
}

LastPatternTree::LastPatternTree ( LastPatternTree &parenttree, unsigned int legindex ) {
  LastLeg &leg = * ( parenttree.legs[legindex] );
    
  LastaddCloseExtensions ( closelegs, parenttree.closelegs, leg.occurrences );
  
  symmetric = parenttree.symmetric;
  // update information used to determine canonical form
  int i, pos;

  int maxleg;
  

  treetuples.reserve ( parenttree.treetuples.size () );
  if ( !parenttree.rootpathstart && (int) legindex < parenttree.secondpathleg ) {
    for ( i = 0; i < (int) parenttree.nextpathstart; i++ )
      treetuples.push_back ( parenttree.treetuples[i] );
    pos = treetuples.size ();
    treetuples.push_back ( leg.tuple );
    nextpathstart = treetuples.size ();
    rootpathstart = 0;

    for ( ; i < (int) parenttree.treetuples.size (); i++ )
      treetuples.push_back ( parenttree.treetuples[i] );

    maxleg = parenttree.secondpathleg;
  }
  else {
    nextpathstart = rootpathstart = parenttree.nextpathstart;
    for ( i = 0; i < (int) parenttree.treetuples.size (); i++ )
      treetuples.push_back ( parenttree.treetuples[i] );
    pos = treetuples.size ();
    treetuples.push_back ( leg.tuple );

    maxleg = parenttree.legs.size ();
  }

  maxdepth = parenttree.maxdepth;

  rootpathrelations.reserve ( leg.tuple.depth + 2 );
  rightmostindexes.reserve ( leg.tuple.depth + 2 );
  if ( parenttree.rootpathstart == 0 && (int) legindex >= parenttree.secondpathleg ) {
    // we're going from the first path to the second
    for ( i = 0; i < (int) leg.tuple.depth; i++ ) {
      rootpathrelations.push_back ( 0 );
      rightmostindexes.push_back ( rootpathstart + i );
    }
    rootpathrelations.push_back ( treetuples[rootpathstart + i].label - leg.tuple.label );
    rightmostindexes.push_back ( treetuples.size () - 1 );
    if ( symmetric && treetuples[maxdepth + 1] == treetuples.back () )
      // in this case, we're copying the previous tree
      nextprefixindex = maxdepth + 2;
    else
      if ( treetuples[rootpathstart + leg.tuple.depth].label == leg.tuple.label )
	nextprefixindex = rootpathstart + leg.tuple.depth + 1;
      else
	nextprefixindex = NONEXTPREFIX;
  }
  else {
    for ( i = 0; i < (int) leg.tuple.depth; i++ ) {
      rootpathrelations.push_back ( parenttree.rootpathrelations[i] );
      rightmostindexes.push_back ( parenttree.rightmostindexes[i] );
    }
    rightmostindexes.push_back ( pos );
    if ( rootpathrelations.size () == 0 || rootpathrelations.back () ==  0 )
      rootpathrelations.push_back ( treetuples[rootpathstart + leg.tuple.depth].label - leg.tuple.label );
    else
      rootpathrelations.push_back ( rootpathrelations.back () );

    if ( parenttree.nextprefixindex != NONEXTPREFIX && parenttree.treetuples[parenttree.nextprefixindex] == leg.tuple )
      nextprefixindex = parenttree.nextprefixindex + 1;
    else {
      if ( leg.tuple.depth < parenttree.rightmostindexes.size () &&
           parenttree.treetuples[parenttree.rightmostindexes[leg.tuple.depth]].label == leg.tuple.label )
        nextprefixindex = parenttree.rightmostindexes[leg.tuple.depth] + 1;
      else
        nextprefixindex = NONEXTPREFIX;
    }
  }

  // ADDED
  fm::last_graphstate->treetuples = &treetuples;
  fm::last_graphstate->closetuples = NULL;
  fm::last_graphstate->startsecondpath = nextpathstart;
    
  if ( nextprefixindex == nextpathstart && symmetric == 1 ) {
    secondpathleg = 0; // THE BUG
    bbrc_extend ( leg.occurrences, MAXEDGELABEL, (unsigned char) NONODE );
    LastaddCloseExtensions ( closelegs, leg.occurrences.number );
    return;
  }

  // determine legs that can be added at leg.tuple.depth + 1 (if possible) (type 1)

  addExtensionLastLegs ( leg.tuple, leg.occurrences );

  int index = legindex;

  if ( nextprefixindex != NONEXTPREFIX && treetuples[nextprefixindex].depth <= leg.tuple.depth ) {
    LastDepth nextprefixdepth = treetuples[nextprefixindex].depth;
    LastEdgeLabel nextprefixlabel = treetuples[nextprefixindex].label;
    LastEdgeLabel lowestlabel = treetuples[rootpathstart + nextprefixdepth].label;
    while ( index < maxleg  &&
            parenttree.legs[index]->tuple.depth > nextprefixdepth )
      index++;
    if ( index < maxleg && parenttree.legs[index]->tuple.depth == nextprefixdepth &&
         nextprefixlabel != lowestlabel ) {
      if ( parenttree.legs[index]->tuple.label == lowestlabel )
        index++;
      while ( index < maxleg &&
              parenttree.legs[index]->tuple.depth == nextprefixdepth &&
	      parenttree.legs[index]->tuple.label < nextprefixlabel )
	index++;
    }
  }

  if ( index == (int) legindex ) {
    LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences );
    if ( legoccurrencesptr )
      addLastLeg ( leg.tuple.connectingnode, leg.tuple.depth, leg.tuple.label, *legoccurrencesptr );
    index++;
  }

  if ( rootpathstart == 0 ) {
    secondpathleg = legs.size (); // THE BUG
    while ( index < (int) parenttree.legs.size () ) {
      if ( index == parenttree.secondpathleg )
        secondpathleg = legs.size ();
      LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, parenttree.legs[index]->tuple.connectingnode, parenttree.legs[index]->occurrences );
      if ( legoccurrencesptr )
        addLastLeg ( parenttree.legs[index]->tuple.connectingnode, parenttree.legs[index]->tuple.depth, parenttree.legs[index]->tuple.label, *legoccurrencesptr );
      index++;
    }
    if ( index == parenttree.secondpathleg )
      secondpathleg = legs.size ();
  }
  else {
    while ( index < (int) parenttree.legs.size () ) {
      LastLegOccurrencesPtr legoccurrencesptr = bbrc_join ( leg.occurrences, parenttree.legs[index]->tuple.connectingnode,  parenttree.legs[index]->occurrences );
      if ( legoccurrencesptr )
        addLastLeg ( parenttree.legs[index]->tuple.connectingnode, parenttree.legs[index]->tuple.depth, parenttree.legs[index]->tuple.label, *legoccurrencesptr );
      index++;
    }
    secondpathleg = legs.size ();
  }
}

GSWalk* LastPatternTree::expand (pair<float, string> max, const int parent_size) {

  assert(parent_size>0);

  fm::last_statistics->patternsize++;
  if ( fm::last_statistics->patternsize > (int) fm::last_statistics->frequenttreenumbers.size () ) {
    fm::last_statistics->frequenttreenumbers.resize ( fm::last_statistics->patternsize, 0 );
    fm::last_statistics->frequentpathnumbers.resize ( fm::last_statistics->patternsize, 0 );
    fm::last_statistics->frequentgraphnumbers.resize ( fm::last_statistics->patternsize, 0 );
  }
  ++fm::last_statistics->frequenttreenumbers[fm::last_statistics->patternsize-1];
  if ( fm::last_statistics->patternsize == ((1<<(sizeof(LastNodeId)*8))-1) ) {
    fm::last_statistics->patternsize--;
    return NULL;
  }
   

  // new siblingwalk
  GSWalk* siblingwalk = new GSWalk();

  // needed for topdown and sibling merge
  vector<int> core_ids; 
  for (int j=0; j<parent_size; j++) core_ids.push_back(j);

  for ( int i=legs.size()-1; i>=0; i-- ) {


    // new current pattern
    GSWalk* gsw = new GSWalk();
    GSWalk* topdown = NULL;

    bool nsign=1;

    float cur_chisq;
    if (fm::last_chisq->active) {
        if (!fm::last_regression) { fm::last_chisq->Calc(legs[i]->occurrences.elements); cur_chisq=fm::last_chisq->p; }
        else                 {    fm::last_ks->Calc(legs[i]->occurrences.elements); cur_chisq=   fm::last_ks->p; }
    }
    

    fm::last_graphstate->insertNode ( legs[i]->tuple.connectingnode, legs[i]->tuple.label, legs[i]->occurrences.maxdegree );
    #ifdef DEBUG
    fm::last_graphstate->print(legs[i]->occurrences.frequency);
    #endif

    #ifdef DEBUG
    fm::last_gsp_out=false;
    string s = fm::last_graphstate->to_s(legs[i]->occurrences.frequency);
    bool diehard=0;
    //if (s.find("N-C-C(-O-C-N)(=C-C)")!=string::npos) { fm::last_die=1; diehard=1; }
    #endif

    if (fm::last_chisq->active) { 
        map<float, map<LastTid, int> > weightmap; 
        for (map<float, LastTid>::iterator it=fm::last_chisq->nr_acts.begin(); it!=fm::last_chisq->nr_acts.end(); it++) {
          weightmap[it->first].insert(make_pair(it->second, 1)); 
        }
//        map<LastTid, int> weightmap_a; each_it(fm::last_chisq->fa_set, set<LastTid>::iterator) { weightmap_a.insert(make_pair((*it),1)); }
//        map<LastTid, int> weightmap_i; each_it(fm::last_chisq->fi_set, set<LastTid>::iterator) { weightmap_i.insert(make_pair((*it),1)); }
        fm::last_graphstate->print(gsw, weightmap); // print to graphstate walk

        if (!fm::last_regression) {
            gsw->activating=fm::last_chisq->activating;
            if (cur_chisq >= fm::last_chisq->sig) {
                nsign=0;
            }
        }
        else {
            gsw->activating=fm::last_ks->activating;
            if (cur_chisq >= fm::last_ks->sig) {
                nsign=0;
            }
        }

    }
    const int gsw_size = gsw->nodewalk.size();

    // !STOP: MERGE TO SIBLINGWALK
    if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Already nodes marked as available 5.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl;exit(1); }
    if (nsign || gsw->activating!=siblingwalk->activating || siblingwalk->hops > fm::last_max_hops) { // empty sw needs no checks
          if (siblingwalk->hops>0) {
              if (siblingwalk->hops>1) {
                    siblingwalk->svd();
              }
              if (fm::last_do_output) {
                  if (!fm::last_console_out) { 
                      ostringstream strstrm;
                      strstrm << siblingwalk;
                      (*fm::last_result) << strstrm.str();
                  }
                  else cout << siblingwalk;
              }
          }
          delete siblingwalk;
          siblingwalk = new GSWalk();
    }
    if (!nsign && ((gsw->activating==siblingwalk->activating) || !siblingwalk->edgewalk.size())) {
        #ifdef DEBUG
        if (fm::last_die) cout << "CR gsw" << endl;
        #endif
        int res=gsw->conflict_resolution(core_ids, siblingwalk);
    }

    if (gsw->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr<<"Error! Still nodes marked as available 5.1. "<<gsw->to_nodes_ex.size()<<" "<<siblingwalk->to_nodes_ex.size()<<endl; exit(1); }

    
    // RECURSE
    if ( ( !fm::last_do_pruning ||  (  fm::last_chisq->u >= fm::last_chisq->sig) ) &&
         (  fm::last_refine_singles || (legs[i]->occurrences.frequency>1) )
       ) {
        LastPatternTree p ( *this, i );
        if (cur_chisq > max.first) { fm::last_updated = true; topdown = p.expand (pair<float, string>(cur_chisq,fm::last_graphstate->to_s(legs[i]->occurrences.frequency)), gsw_size); }
        else topdown = p.expand (max, gsw_size);
    }

    // merge to siblingwalk
    if (topdown != NULL) {
       if (topdown->edgewalk.size()) {

            #ifdef DEBUG
            if (fm::last_die) {
                cout << "TOPDOWN2 BEGIN " << core_ids.size() << endl;
                cout << topdown ;
                cout << "--result--" << endl;
                cout << siblingwalk ;
            }
            #endif

            if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Already nodes marked as available 5.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }
            // STOP: OUTPUT TOPDOWN
            if (nsign || siblingwalk->activating!=topdown->activating) { 
                #ifdef DEBUG
                if (fm::last_die) cout << "STOP CRITERIUM at CHI " << cur_chisq << endl;
                #endif
                if (topdown->hops>0) { 
                    if (topdown->hops>1) { 
                        topdown->svd();
                    }
                    if (fm::last_do_output) {
                        if (!fm::last_console_out) { 
                            ostringstream strstrm;
                            strstrm << topdown;
                            (*fm::last_result) << strstrm.str();
                        }
                        else cout << topdown;
                    }
                }
            }
            // ELSE: MERGE TO SIBLINGWALK
            else {
                int res=topdown->conflict_resolution(core_ids, siblingwalk); 
            }
            if (topdown->to_nodes_ex.size() || siblingwalk->to_nodes_ex.size()) { cerr << "Error! Still nodes marked as available 5.2. " << topdown->to_nodes_ex.size() << " " << siblingwalk->to_nodes_ex.size() <<  endl; exit(1); }

            #ifdef DEBUG
            if (fm::last_die) {
                cout << "TOPDOWN2 END " << core_ids.size() << endl;
                cout << topdown ;
                cout << "--result--" << endl;
                cout << siblingwalk ;
            }
            #endif

       }
    }
    
    fm::last_graphstate->deleteNode ();
    delete topdown;
    delete gsw;
    #ifdef DEBUG
    if (diehard==1) { 
       cerr << "DYING HARD!" << endl;
       exit(0);
    }
    #endif

  }

  #ifdef DEBUG  
  if (!legs.size()) cout << fm::last_graphstate->sep() << endl;
  #endif
  
  fm::last_statistics->patternsize--;
  return siblingwalk;

}



LastPatternTree::~LastPatternTree () {
  for ( int i = 0; i < (int) legs.size (); i++ )
    delete legs[i];
  for ( int i = 0; i < (int) closelegs.size (); i++ )
    delete closelegs[i];
}

/*
ostream &operator<< ( ostream &stream, LastTuple &tuple ) {
  LastDatabaseLastEdgeLabel edgelabel = database->edgelabels[fm::last_database->edgelabelsindexes[tuple.label]];
  stream << "(" << tuple.depth << ","
         << fm::last_database->nodelabels[edgelabel.fromnodelabel].inputlabel << "-"
         << edgelabel.inputedgelabel << "-"
         << fm::last_database->nodelabels[edgelabel.tonodelabel].inputlabel << "[" << (int) tuple.label << "])";

  return stream;
}
*/

void LastPatternTree::checkIfIndeedNormal () {
  int i = 0, j = nextpathstart;
  bool equal = true;
  for ( ; i <= (int) maxdepth && equal; i++, j++ )
    equal = ( treetuples[i] == treetuples[j] );
  if ( !equal && treetuples[i-1].label > treetuples[j-1].label ) {
    cout << "NOT NORMAL: first path is higher than second" << endl;
  }
  for ( i = maxdepth + 1; i < (int) nextpathstart; i++ ) {
    if ( treetuples[i].depth == maxdepth &&
         treetuples[i].label < treetuples[treetuples[i].depth].label )
      cout << "NOT NORMAL: lower path than possible (1)" << endl;
  }
  for ( i = nextpathstart + maxdepth + 1; i < (int) treetuples.size (); i++ ) {
    if ( treetuples[i].depth == 0 ) {
      j = i; i = nextpathstart;
      equal = true;
      for ( ; j < (int) treetuples.size () && equal; i++, j++ )
        equal = ( treetuples[i] == treetuples[j] );
      if ( !equal && treetuples[j-1].depth == maxdepth && treetuples[i-1].label > treetuples[j-1].label ) {
        cout << "NOT NORMAL: third path is lower than the second" << endl;
      }
      break;
    }
    if ( treetuples[i].depth == maxdepth && 
         treetuples[i].label < treetuples[treetuples[i].depth + nextpathstart].label )
      cout << "NOT NORMAL: lower path than possible (2)" << endl;
  }
}

ostream& operator<< (ostream& os, GSWalk* gsw) {
    static int gsw_counter=0;

    #ifndef DEBUG
    if (gsw->edgewalk.size()) {
        gsw_counter++;
        os << "    <graph id=\"" << gsw_counter << "\" edgedefault=\"undirected\">" << endl;
        os << "        <data key=\"act\">" << gsw->activating << "</data>" << endl;
        os << "        <data key=\"hops\">" << gsw->hops << "</data>" << endl;
    }

    for(vector<GSWNode>::iterator it=gsw->nodewalk.begin(); it!=gsw->nodewalk.end(); it++) {
        os << "        <node id=\"" << distance(gsw->nodewalk.begin(), it) << "\">" << endl;
        string labels;
        for (set<InputLastNodeLabel>::iterator it2=it->labs.begin(); it2!=it->labs.end(); it2++) {
            if (it2!=it->labs.begin()) labels.append(" ");
            labels.append(to_string(*it2));
        }
        os << "            <data key=\"lab_n\">" << labels << "</data>" << endl;
        os << "        </node>" << endl;
    }

    for (map<int, map<int, GSWEdge> >::iterator it=gsw->edgewalk.begin(); it!=gsw->edgewalk.end(); it++) {

        for(map<int,GSWEdge>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            os << "        <edge source=\"" << it->first << "\" target=\"" << it2->first << "\">" << endl;

            // from and to
            string labels;
            for (set<InputLastEdgeLabel>::iterator it3=it2->second.labs.begin(); it3!=it2->second.labs.end(); it3++) {
                if (it3!=it2->second.labs.begin()) labels.append(" ");
                labels.append(to_string(*it3));
            }
            os << "            <data key=\"lab_e\">" << labels << "</data>" << endl;
            os << "            <data key=\"weight\">" << it2->second.discrete_weight << "</data>" << endl;
            os << "            <data key=\"del\">" << it2->second.deleted << "</data>" << endl;
            os << "        </edge>" << endl;
        }
    }

    if (gsw->edgewalk.size()) {
        os << "    </graph>" << endl;
        os << endl;
    }
    #endif
    
    #ifdef DEBUG
    if (gsw->edgewalk.size()) {
        gsw_counter++;
        //os << "t # " << gsw_counter << " " << gsw->activating << " " << gsw->hops << endl;
    }

    for(vector<GSWNode>::iterator it=gsw->nodewalk.begin(); it!=gsw->nodewalk.end(); it++) {
        os << distance(gsw->nodewalk.begin(), it);
        os << " < ";
        for (set<InputLastNodeLabel>::iterator it2=it->labs.begin(); it2!=it->labs.end(); it2++) {
            os << *it2 << " ";
        }
        os << ">";
        os << endl;
    }

    for (map<int, map<int, GSWEdge> >::iterator it=gsw->edgewalk.begin(); it!=gsw->edgewalk.end(); it++) {
        if (!it->second.size()) os << it->first << endl;

        for(map<int,GSWEdge>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {

            // from and to
            os << it->first << " " << it2->first; 

            os << " < ";
            for (set<InputLastEdgeLabel>::iterator it3=it2->second.labs.begin(); it3!=it2->second.labs.end(); it3++) {
                os << *it3 << " ";
            }
            os << ">";

            /*
            int count=0;
            for (map<LastTid, int>::iterator it3=it2->second.a.begin(); it3!=it2->second.a.end(); it3++) {
                count = count + it3->second;
            }
            for (map<LastTid, int>::iterator it3=it2->second.i.begin(); it3!=it2->second.i.end(); it3++) {
                count = count + it3->second;
            }
            os << " " << count;
            */
            os << " " << it2->second.discrete_weight;

            if (it2->second.deleted) os << " [";
            if (it2->second.deleted) os << "D";
            if (it2->second.deleted) os << "] ";

            os << endl;
        }
    }

    if (gsw->edgewalk.size()) {
        os << endl;
    }
    #endif

    return os;
};


ostream& operator<< (ostream& os, GSWEdge* gswe) {
    typedef map<LastTid,int> mmap;
    typedef map<float, map<LastTid,int> >mmmap;
    os << "To: " << gswe->to;
    os << " Labs: <";
    each_it(gswe->labs, set<InputLastEdgeLabel>::iterator) {
        os << *it << " ";
    }
    os << "> ";
    for (mmmap::iterator it=gswe->m.begin(); it!=gswe->m.end(); it++) {
      os << "> ";
      os << "<";
      for (mmap::iterator it2=it->second.begin(); it2!=it->second.end(); it2++) {
        os << it2->first << "->" << it2->second << " ";
      }
    }
    os << "> ";
}

ostream& operator<< (ostream& os, GSWNode* gswn) {
    typedef map<LastTid,int> mmap;
    os << " Labs: <";
    each_it(gswn->labs, set<InputLastEdgeLabel>::iterator) {
        os << *it << " ";
    }
    os << "> ";
}
