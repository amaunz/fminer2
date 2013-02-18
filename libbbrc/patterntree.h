// patterntree.h
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

#ifndef PATTERNTREE_H
#define PATTERNTREE_H
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

#include "misc.h"
#include "database.h"
#include "legoccurrence.h"
#include "path.h"
#include "closeleg.h"

using namespace std;

extern bool updated;
extern float minfreq;

struct BbrcTuple {
  BbrcDepth depth;
  BbrcEdgeLabel label;
  BbrcNodeId connectingnode;
  
  BbrcTuple ( BbrcDepth depth, BbrcEdgeLabel label, BbrcNodeId connectingnode ): 
    depth ( depth ), label ( label ), connectingnode ( connectingnode ) { }
  BbrcTuple () { }

  friend bool operator< ( BbrcTuple &a, BbrcTuple &b ) { return a.depth > b.depth || ( a.depth == b.depth && a.label < b.label ); }
  friend bool operator== ( BbrcTuple &a, BbrcTuple &b ) { return a.depth == b.depth && a.label == b.label; }
  friend ostream &operator<< ( ostream &stream, BbrcTuple &tuple );
};

struct BbrcLeg {
  BbrcTuple tuple;
  BbrcLegOccurrences occurrences;
};

typedef BbrcLeg *BbrcLegPtr;

class BbrcPatternTree {
  public:
    BbrcPatternTree ( BbrcPath &path, unsigned int legindex );
    ~BbrcPatternTree ();
    void expand (pair<float, string> max);
    vector<BbrcLegPtr> legs; // pointers used to avoid copy-constructor during a resize of the vector
  private:
    void checkIfIndeedNormal ();
    /* inline */ void addExtensionBbrcLegs ( BbrcTuple &tuple, BbrcLegOccurrences &legoccurrences );
    /* inline */ void addBbrcLeg ( const BbrcNodeId connectingnode, const int depth, const BbrcEdgeLabel edgelabel, BbrcLegOccurrences &legoccurrences );
    /* inline */ void addLeftBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, int &i, BbrcDepth olddepth, BbrcEdgeLabel lowestlabel, int leftend, int edgesize2 );
    /* inline */ int addLeftBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, BbrcTuple &tuple, unsigned int legindex, int leftend, int edgesize2 );
    /* inline */ void addRightBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, int &i, BbrcDepth olddepth, BbrcEdgeLabel lowestlabel, int rightstart, int nodesize2 );
    /* inline */ int addRightBbrcLegs ( BbrcPath &path, BbrcPathBbrcLeg &leg, BbrcTuple &tuple, unsigned int legindex, int rightstart, int nodesize2 );
    BbrcPatternTree ( BbrcPatternTree &parenttree, unsigned int legindex );
    vector<BbrcTuple> treetuples;
    vector<BbrcNodeId> rightmostindexes;
    vector<short> rootpathrelations;
    unsigned int nextprefixindex;
    unsigned int rootpathstart;
    unsigned int nextpathstart;
    unsigned int maxdepth;
    int symmetric; // 0 == not symmetric, 1 == symmetric, even length path, 2 == symmetric, odd length path
    int secondpathleg;
    vector<BbrcCloseBbrcLegPtr> closelegs;
    friend ostream &operator<< ( ostream &stream, BbrcPatternTree &patterntree );
#ifdef GRAPH_OUTPUT
    friend void fillMatrix ( int **A, int &nextnode, int rootnode, BbrcNodeLabel rootlabel, 
                  int startpos, int endpos, BbrcPatternTree &patterntree );
    BbrcNodeLabel tree1rootlabel, tree2rootlabel;
    BbrcEdgeLabel rootpathlabel;
#endif
};

#define NONEXTPREFIX ((unsigned int) -1)

#endif
