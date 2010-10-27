// patterntree.cpp
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

#include "patterntree.h"
#include "graphstate.h"

namespace fm {
    extern unsigned int minfreq;
    extern bool do_backbone;
    extern bool updated;
    extern bool adjust_ub;
    extern bool do_pruning;
    extern bool console_out;
    extern bool refine_singles;
    extern bool do_output;
    extern bool bbrc_sep;
    extern bool regression;

    extern BbrcDatabase* database;
    extern ChisqConstraint* chisq;
    extern KSConstraint* ks;
    extern vector<string>* result;
    extern BbrcStatistics* statistics;
    extern BbrcGraphState* graphstate;
    extern BbrcBbrcLegOccurrences* legoccurrences;

    extern vector<BbrcBbrcLegOccurrences> Bbrccandidatelegsoccurrences; 
}

int maxsize = ( 1 << ( sizeof(BbrcNodeId)*8 ) ) - 1; // safe default for the largest allowed pattern

inline void BbrcPatternTree::addBbrcLeg ( BbrcNodeId connectingnode, const int depth, const BbrcEdgeLabel edgelabel, BbrcBbrcLegOccurrences &legoccurrences ) {
  BbrcLegPtr leg = new BbrcLeg;
  leg->tuple.depth = depth;
  leg->tuple.label = edgelabel;
  leg->tuple.connectingnode = connectingnode;
  store ( leg->occurrences, legoccurrences );
  legs.push_back ( leg );
}

// this function assumes that the extension tuple is already added at the back of the queue,
// and the equivalency information has been filled in.
void BbrcPatternTree::addExtensionBbrcLegs ( BbrcTuple &tuple, BbrcBbrcLegOccurrences &legoccurrences ) {
  if ( legoccurrences.maxdegree == 1 )
    return;
  if ( tuple.depth == maxdepth ) {
    extend ( legoccurrences, MAXEDGELABEL, (unsigned char) NONODE );
    BbrcaddCloseExtensions ( closelegs, legoccurrences.number );
    return;
  }
  BbrcEdgeLabel minlabel = NOEDGELABEL, neglect = '\0', pathlowestlabel = treetuples[tuple.depth + 1 + rootpathstart].label;

  if ( nextprefixindex != NONEXTPREFIX ) {
    if ( treetuples[nextprefixindex].depth <= tuple.depth ) {
      // heuristic saving
      extend ( legoccurrences, MAXEDGELABEL, (unsigned char) NONODE );
      BbrcaddCloseExtensions ( closelegs, legoccurrences.number );
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
      extend ( legoccurrences, MAXEDGELABEL, (unsigned char) NONODE );
      BbrcaddCloseExtensions ( closelegs, legoccurrences.number );
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
    extend ( legoccurrences, minlabel, neglect );
  else
    extend ( legoccurrences );

  if ( fm::Bbrccandidatelegsoccurrences[pathlowestlabel].frequency >= fm::minfreq )
    // this is the first possible extension, as we force this label to be the lowest!
    addBbrcLeg ( fm::graphstate->lastNode (), tuple.depth + 1, pathlowestlabel, fm::Bbrccandidatelegsoccurrences[pathlowestlabel] );

  for ( int i = 0; (unsigned) i < fm::Bbrccandidatelegsoccurrences.size (); i++ ) {
    if ( fm::Bbrccandidatelegsoccurrences[i].frequency >= fm::minfreq && i != pathlowestlabel )
      addBbrcLeg ( fm::graphstate->lastNode (), tuple.depth + 1, i, fm::Bbrccandidatelegsoccurrences[i] );
  }

  BbrcaddCloseExtensions ( closelegs, legoccurrences.number );
}

void BbrcPatternTree::addLeftBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, int &i, BbrcDepth olddepth, BbrcEdgeLabel lowestlabel, int leftend, int edgesize2 ) {
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
          BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, path.legs[i2]->tuple.connectingnode, path.legs[i2]->occurrences );
          if ( legoccurrencesptr )
            addBbrcLeg ( path.legs[i2]->tuple.connectingnode, edgesize2 - olddepth, path.legs[i2]->tuple.edgelabel, *legoccurrencesptr );
          break;
        }
        i2++;
      }
    }
    // skip lowest label tuples, as they have already been moved to the front...
    if ( path.legs[i]->tuple.edgelabel != lowestlabel ) {
      BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, path.legs[i]->tuple.connectingnode, path.legs[i]->occurrences );
      if ( legoccurrencesptr )
        addBbrcLeg ( path.legs[i]->tuple.connectingnode, edgesize2 - olddepth, path.legs[i]->tuple.edgelabel, *legoccurrencesptr );
    }
  }
}

int BbrcPatternTree::addLeftBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, BbrcTuple &tuple, unsigned int legindex, int leftend, int edgesize2 ) {
  int i;

  BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences );
  if ( legoccurrencesptr )
    addBbrcLeg ( leg.tuple.connectingnode, tuple.depth, tuple.label, *legoccurrencesptr );

  // the easy part - the extensions of the left side of the original path

  if ( rootpathrelations.back () == 0 && tuple.depth != maxdepth ) {
    // we force the label on the path to be the lowest label, so all labels which are actually
    // lower must be added as they are higher here, deepest nodes are an exception;
    // no path may be lower than the current path

    for ( i = legindex - 1; i >= 0 && path.legs[i]->tuple.depth == leg.tuple.depth; i-- );
    BbrcEdgeLabel lowestlabel = path.edgelabels[leg.tuple.depth - 1];
    for ( i++; i < (int) legindex; i++ ) {
      if ( path.legs[i]->tuple.edgelabel != lowestlabel ) {
        legoccurrencesptr = join ( leg.occurrences, path.legs[i]->tuple.connectingnode, path.legs[i]->occurrences );
        if ( legoccurrencesptr )
          addBbrcLeg ( path.legs[i]->tuple.connectingnode, tuple.depth, path.legs[i]->tuple.edgelabel, *legoccurrencesptr );
      }
    }
  }

  // the order of the other extensions is almost correct - except that we must move
  // extensions at depth d which are for the label of the path at depth d to the front of
  // all legs of depth d...

  i = legindex + 1;
  addLeftBbrcLegs ( path, leg, i, leg.tuple.depth, path.edgelabels[leg.tuple.depth - 1], leftend, edgesize2 );

  return i;
}

void BbrcPatternTree::addRightBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, int &i, BbrcDepth olddepth, BbrcEdgeLabel lowestlabel, int rightstart, int nodesize2 ) {
  int i2 = i + 1, k;
  BbrcBbrcLegOccurrencesPtr legoccurrencesptr;

  while ( i >= 0 && (int) path.legs[i]->tuple.depth >= rightstart ) {
    // we encounter a new depth, do all the extensions of the previous depths...
    if ( path.legs[i]->tuple.depth != olddepth ) {
      for ( k = i + 1; k < i2; k++ ) {
        if ( path.legs[k]->tuple.edgelabel != lowestlabel ) {
          legoccurrencesptr = join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
          if ( legoccurrencesptr )
            addBbrcLeg ( path.legs[k]->tuple.connectingnode, path.legs[k]->tuple.depth - nodesize2, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
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
      legoccurrencesptr = join ( leg.occurrences, path.legs[i]->tuple.connectingnode, path.legs[i]->occurrences );
      if ( legoccurrencesptr )
        addBbrcLeg ( path.legs[i]->tuple.connectingnode, path.legs[i]->tuple.depth - nodesize2, path.legs[i]->tuple.edgelabel, *legoccurrencesptr );
    }
    i--;
  }
  // some tuples may not have been checked yet
  for ( k = i + 1; k < i2; k++ ) {
    if ( path.legs[k]->tuple.edgelabel != lowestlabel ) {
      legoccurrencesptr = join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
      if ( legoccurrencesptr )
        addBbrcLeg ( path.legs[k]->tuple.connectingnode, path.legs[k]->tuple.depth - nodesize2, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
    }
  }
}

int BbrcPatternTree::addRightBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, BbrcTuple &tuple, unsigned int legindex, int rightstart, int nodesize2 ) {
  int i;
  BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences );
  if ( legoccurrencesptr )
    addBbrcLeg ( leg.tuple.connectingnode, tuple.depth, tuple.label, *legoccurrencesptr );

  // other extensions at the right path

  i = legindex - 1;
  while ( i >= 0 && path.legs[i]->tuple.depth == leg.tuple.depth )
    i--;

  // brothers of the new leg
  if ( rootpathrelations.back () == 0 && tuple.depth != maxdepth )
    for ( int j = i + 1; j < (int) legindex; j++ ) {
      BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, path.legs[j]->tuple.connectingnode, path.legs[j]->occurrences );
      if ( legoccurrencesptr )
        addBbrcLeg ( path.legs[j]->tuple.connectingnode, tuple.depth, path.legs[j]->tuple.edgelabel, *legoccurrencesptr );
    }
  BbrcEdgeLabel lowestlabel = path.edgelabels[leg.tuple.depth];
  for ( int j = legindex + 1; j < (int) path.legs.size () && path.legs[j]->tuple.depth == leg.tuple.depth; j++ ) {
    if ( path.legs[j]->tuple.edgelabel != lowestlabel ) {
      BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, path.legs[j]->tuple.connectingnode, path.legs[j]->occurrences );
      if ( legoccurrencesptr )
        addBbrcLeg ( path.legs[j]->tuple.connectingnode, tuple.depth, path.legs[j]->tuple.edgelabel, *legoccurrencesptr );
    }
  }

  // candidates at a lower depth
  if ( i >= 0 )
    addRightBbrcLegs ( path, leg, i, leg.tuple.depth, lowestlabel, rightstart, nodesize2 );

  return i;
}

BbrcPatternTree::BbrcPatternTree ( BbrcPath &path, unsigned int legindex ) {
  BbrcPathBbrcLeg &leg = (*path.legs[legindex]);
  
  maxdepth = path.edgelabels.size () / 2 - 1;
  int leftwalk, leftstart, rightwalk, rightstart;
  BbrcBbrcLegOccurrencesPtr legoccurrencesptr;

  BbrcaddCloseExtensions ( closelegs, path.closelegs, leg.occurrences );

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
    fm::graphstate->nasty = ( leftwalk == -1 );
    if ( leftwalk == -1 || path.edgelabels[leftwalk] < path.edgelabels[rightwalk] ) {
      // left part of the path should be the first path in the tree

      treetuples.reserve ( path.edgelabels.size () + 1 );
      int i, j;
      for ( i = leftstart, j = 0; i >= 0; i--, j++ ) {
        vector_push_back ( BbrcTuple, treetuples, tuple );
        tuple.label = path.edgelabels[i];
        tuple.depth = j;
      }

      int leftend = nodesize2 - 1;

      if ( (int) leg.tuple.depth <= leftend ) {
        // extension is an extension of the first path

        // add extension tuple
        vector_push_back ( BbrcTuple, treetuples, tuple );
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
        addExtensionBbrcLegs ( tuple, leg.occurrences );

        // add second path
        for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
          vector_push_back ( BbrcTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }

        // here, we have determined the expansion tuple sequence of the tree, and the normal
        // form information. Now we have to fill in the legs at other positions than below
        // the new node. This is complicated as the occurrence list order of the path
        // has to be changed...

        i = addLeftBbrcLegs ( path, leg, tuple, legindex, leftend, edgesize2 );

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
              legoccurrencesptr = join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
              if ( legoccurrencesptr ) {
                addBbrcLeg ( path.legs[k]->tuple.connectingnode, path.legs[k]->tuple.depth - nodesize2, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
              }
            }
          }
          while ( j >= i && path.legs[j]->tuple.depth == path.nodelabels.size () - 2 )
            j--;
        }

        if ( j >= i )
          addRightBbrcLegs ( path, leg, j, NODEPTH, NOEDGELABEL, rightstart, nodesize2 );

        // we're done!

      }
      else {
        // extension is an extension of the second path

        rootpathstart = nextpathstart = treetuples.size ();

        // add second path
        for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
          vector_push_back ( BbrcTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }
        // add extension tuple
        vector_push_back ( BbrcTuple, treetuples, tuple );
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
        addExtensionBbrcLegs ( tuple, leg.occurrences );

        // add the leg itself

        addRightBbrcLegs ( path, leg, tuple, legindex, rightstart, nodesize2 );

        secondpathleg = legs.size ();
        // we're done
      }
    }
    else {
      // right part of the tree is the first path in the tree

      treetuples.reserve ( path.edgelabels.size () + 1 );
      int i, j;
      for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
        vector_push_back ( BbrcTuple, treetuples, tuple );
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
          vector_push_back ( BbrcTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }
        // add extension tuple
        vector_push_back ( BbrcTuple, treetuples, tuple );
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
        addExtensionBbrcLegs ( tuple, leg.occurrences );

        // add the leg itself, and all legs on the left path above in the tree

        addLeftBbrcLegs ( path, leg, tuple, legindex, leftend, edgesize2 );

        secondpathleg = legs.size ();

        // we're done
      }
      else {
        // the right path is the lowest, and we're adding a leg in that path
        // add extension tuple

        int i, j;

        vector_push_back ( BbrcTuple, treetuples, tuple );
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
        addExtensionBbrcLegs ( tuple, leg.occurrences );

        // add second path
        for ( i = leftstart, j = 0; i >= 0; i--, j++ ) {
          vector_push_back ( BbrcTuple, treetuples, tuple );
          tuple.label = path.edgelabels[i];
          tuple.depth = j;
        }

        // here, we have determined the expansion tuple sequence of the tree, and the normal
        // form information. Now we have to fill in the legs at other positions than below
        // the new node. This is complicated as the occurrence list order of the path
        // has to be changed...


        i = addRightBbrcLegs ( path, leg, tuple, legindex, rightstart, nodesize2 );

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
          addLeftBbrcLegs ( path, leg, j, NODEPTH, NOEDGELABEL, leftend, edgesize2 );

        // we're done!
      }
    }
  }
  else {

    treetuples.reserve ( path.edgelabels.size () + 1 );
    int i, j;
    for ( i = leftstart, j = 0; i >= 0; i--, j++ ) {
      vector_push_back ( BbrcTuple, treetuples, tuple );
      tuple.label = path.edgelabels[i];
      tuple.depth = j;
    }

    // extension is always an extension of the first path, unless the extension
    // is performed on the middle node of a string of odd number of nodes

    // add extension tuple
    vector_push_back ( BbrcTuple, treetuples, tuple );
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
    addExtensionBbrcLegs ( tuple, leg.occurrences );

    // add second path
    for ( i = rightstart, j = 0; i < (int) path.edgelabels.size (); i++, j++ ) {
      vector_push_back ( BbrcTuple, treetuples, tuple );
      tuple.label = path.edgelabels[i];
      tuple.depth = j;
    }

    // here, we have determined the expansion tuple sequence of the tree, and the normal
    // form information. Now we have to fill in the legs at other positions than below
    // the new node. This is complicated as the occurrence list order of the path
    // has to be changed...

    i = addLeftBbrcLegs ( path, leg, tuple, legindex, rightstart - 1, edgesize2 );

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
        BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
        if ( legoccurrencesptr )
          addBbrcLeg ( path.legs[k]->tuple.connectingnode, tuple.depth, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
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
          BbrcEdgeLabel lowestlabel = path.edgelabels[targetdepth];
          while ( j >= i &&
	          (int) path.legs[j]->tuple.depth == targetdepth &&
	  	  path.legs[j]->tuple.edgelabel >= tuple.label )
            j--;
          for ( int k = j + 1; k <= j2; k++ )
            if ( path.legs[k]->tuple.edgelabel != lowestlabel ) {
              BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, path.legs[k]->tuple.connectingnode, path.legs[k]->occurrences );
              if ( legoccurrencesptr )
                addBbrcLeg ( path.legs[k]->tuple.connectingnode, tuple.depth, path.legs[k]->tuple.edgelabel, *legoccurrencesptr );
          }
          while ( j >= i && (int) path.legs[j]->tuple.depth == targetdepth )
            j--;
        }
      }
    }

    if ( j >= i )
      addRightBbrcLegs ( path, leg, j, NODEPTH, NOEDGELABEL, rightstart, nodesize2 );
        // we're done!

    // symmetry
  }
  
  // ADDED
  fm::graphstate->backbonelength = path.nodelabels.size ();
  if ( fm::graphstate->backbonelength % 2 == 0 )
    fm::graphstate->bicenterlabel = path.edgelabels [ fm::graphstate->backbonelength / 2 - 1 ];
  else
    fm::graphstate->centerlabel = path.nodelabels [ ( fm::graphstate->backbonelength - 1 ) / 2 ];
  fm::graphstate->nasty = false;
  fm::graphstate->treetuples = &treetuples;
  fm::graphstate->closetuples = NULL;
  fm::graphstate->startsecondpath = nextpathstart;
}

BbrcPatternTree::BbrcPatternTree ( BbrcPatternTree &parenttree, unsigned int legindex ) {
  BbrcLeg &leg = * ( parenttree.legs[legindex] );
    
  BbrcaddCloseExtensions ( closelegs, parenttree.closelegs, leg.occurrences );
  
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
  fm::graphstate->treetuples = &treetuples;
  fm::graphstate->closetuples = NULL;
  fm::graphstate->startsecondpath = nextpathstart;
    
  if ( nextprefixindex == nextpathstart && symmetric == 1 ) {
    secondpathleg = 0; // THE BUG
    extend ( leg.occurrences, MAXEDGELABEL, (unsigned char) NONODE );
    BbrcaddCloseExtensions ( closelegs, leg.occurrences.number );
    return;
  }

  // determine legs that can be added at leg.tuple.depth + 1 (if possible) (type 1)

  addExtensionBbrcLegs ( leg.tuple, leg.occurrences );

  int index = legindex;

  if ( nextprefixindex != NONEXTPREFIX && treetuples[nextprefixindex].depth <= leg.tuple.depth ) {
    BbrcDepth nextprefixdepth = treetuples[nextprefixindex].depth;
    BbrcEdgeLabel nextprefixlabel = treetuples[nextprefixindex].label;
    BbrcEdgeLabel lowestlabel = treetuples[rootpathstart + nextprefixdepth].label;
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
    BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences );
    if ( legoccurrencesptr )
      addBbrcLeg ( leg.tuple.connectingnode, leg.tuple.depth, leg.tuple.label, *legoccurrencesptr );
    index++;
  }

  if ( rootpathstart == 0 ) {
    secondpathleg = legs.size (); // THE BUG
    while ( index < (int) parenttree.legs.size () ) {
      if ( index == parenttree.secondpathleg )
        secondpathleg = legs.size ();
      BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, parenttree.legs[index]->tuple.connectingnode, parenttree.legs[index]->occurrences );
      if ( legoccurrencesptr )
        addBbrcLeg ( parenttree.legs[index]->tuple.connectingnode, parenttree.legs[index]->tuple.depth, parenttree.legs[index]->tuple.label, *legoccurrencesptr );
      index++;
    }
    if ( index == parenttree.secondpathleg )
      secondpathleg = legs.size ();
  }
  else {
    while ( index < (int) parenttree.legs.size () ) {
      BbrcBbrcLegOccurrencesPtr legoccurrencesptr = join ( leg.occurrences, parenttree.legs[index]->tuple.connectingnode,  parenttree.legs[index]->occurrences );
      if ( legoccurrencesptr )
        addBbrcLeg ( parenttree.legs[index]->tuple.connectingnode, parenttree.legs[index]->tuple.depth, parenttree.legs[index]->tuple.label, *legoccurrencesptr );
      index++;
    }
    secondpathleg = legs.size ();
  }
}

void BbrcPatternTree::expand (pair<float, string> max) {
  fm::statistics->patternsize++;
  if ( fm::statistics->patternsize > (int) fm::statistics->frequenttreenumbers.size () ) {
    fm::statistics->frequenttreenumbers.resize ( fm::statistics->patternsize, 0 );
    fm::statistics->frequentpathnumbers.resize ( fm::statistics->patternsize, 0 );
    fm::statistics->frequentgraphnumbers.resize ( fm::statistics->patternsize, 0 );
  }
  ++fm::statistics->frequenttreenumbers[fm::statistics->patternsize-1];
  if ( fm::statistics->patternsize == ((1<<(sizeof(BbrcNodeId)*8))-1) ) {
    fm::statistics->patternsize--;
    return;
  }
    
  if (fm::do_backbone && (legs.size()==0)) {
    if (fm::updated)
        if (fm::do_output) {
            if (!fm::console_out) { 
               (*fm::result) << max.second;
            }
            else {
                cout << max.second;
            }
        }
        fm::updated = false;
  }

  
 
  for ( int i = legs.size()-1; i >= 0; i-- ) {

    // Calculate chisq
    if (fm::chisq->active) { 
        if (!fm::regression) fm::chisq->Calc(legs[i]->occurrences.elements);
        else fm::ks->Calc(legs[i]->occurrences.elements);
    }

    // GRAPHSTATE
    fm::graphstate->insertNode ( legs[i]->tuple.connectingnode, legs[i]->tuple.label, legs[i]->occurrences.maxdegree );

    // immediate output for all patterns
    if (fm::do_output && !fm::do_backbone) {
       if (!fm::console_out) (*fm::result) << fm::graphstate->to_s(legs[i]->occurrences.frequency);
       else fm::graphstate->print(legs[i]->occurrences.frequency);
    }

    // RECURSE
    float cmax = maxi ( maxi ( fm::chisq->sig, max.first ), fm::chisq->p );

    if ( ( !fm::do_pruning || 
               (
                 (  !fm::adjust_ub && (fm::chisq->u >= fm::chisq->sig) ) || 
                 (   fm::adjust_ub && (fm::chisq->u >= cmax) )
               )
             ) &&
         (
            fm::refine_singles || (legs[i]->occurrences.frequency>1)
         )
    
    ) {   // UB-PRUNING

        BbrcPatternTree p ( *this, i );

        if (!fm::regression) {
            if (fm::chisq->p > max.first) { fm::updated = true; p.expand (pair<float, string>(fm::chisq->p,fm::graphstate->to_s(legs[i]->occurrences.frequency))); }
            else p.expand (max);
        }
        else {
            if (fm::ks->p > max.first) { fm::updated = true; p.expand (pair<float, string>(fm::ks->p,fm::graphstate->to_s(legs[i]->occurrences.frequency))); }
            else p.expand (max);
        }
    }
    else {
        if (fm::do_backbone && fm::updated) {
            if (fm::do_output) {
                if (!fm::console_out) {
                    (*fm::result) << max.second;
                }
                else {
                    cout << max.second;
                }
            }
            fm::updated = false;
        }
    }

    fm::graphstate->deleteNode ();

  }

  if (fm::bbrc_sep && !fm::do_backbone && (legs.size()==0)) {
      if (fm::do_output) {
          if (!fm::console_out && fm::result->size() && (fm::result->back()!=fm::graphstate->sep())) (*fm::result) << fm::graphstate->sep();
          //else cout << fm::graphstate->sep() << endl;
      }
  }

  fm::statistics->patternsize--;

}



BbrcPatternTree::~BbrcPatternTree () {
  for ( int i = 0; i < (int) legs.size (); i++ )
    delete legs[i];
  for ( int i = 0; i < (int) closelegs.size (); i++ )
    delete closelegs[i];
}

/*
ostream &operator<< ( ostream &stream, BbrcTuple &tuple ) {
  BbrcDatabaseBbrcEdgeLabel edgelabel = database->edgelabels[fm::database->edgelabelsindexes[tuple.label]];
  stream << "(" << tuple.depth << ","
         << fm::database->nodelabels[edgelabel.fromnodelabel].inputlabel << "-"
         << edgelabel.inputedgelabel << "-"
         << fm::database->nodelabels[edgelabel.tonodelabel].inputlabel << "[" << (int) tuple.label << "])";

  return stream;
}
*/

void BbrcPatternTree::checkIfIndeedNormal () {
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
