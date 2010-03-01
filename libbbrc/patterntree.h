// patterntree.h
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
extern unsigned int minfreq;

struct Tuple {
  Depth depth;
  EdgeLabel label;
  NodeId connectingnode;
  
  Tuple ( Depth depth, EdgeLabel label, NodeId connectingnode ): 
    depth ( depth ), label ( label ), connectingnode ( connectingnode ) { }
  Tuple () { }

  friend bool operator< ( Tuple &a, Tuple &b ) { return a.depth > b.depth || ( a.depth == b.depth && a.label < b.label ); }
  friend bool operator== ( Tuple &a, Tuple &b ) { return a.depth == b.depth && a.label == b.label; }
  friend ostream &operator<< ( ostream &stream, Tuple &tuple );
};

struct Leg {
  Tuple tuple;
  LegOccurrences occurrences;
};

typedef Leg *LegPtr;

class PatternTree {
  public:
    PatternTree ( Path &path, unsigned int legindex );
    ~PatternTree ();
    void expand (pair<float, string> max);
    vector<LegPtr> legs; // pointers used to avoid copy-constructor during a resize of the vector
  private:
    void checkIfIndeedNormal ();
    /* inline */ void addExtensionLegs ( Tuple &tuple, LegOccurrences &legoccurrences );
    /* inline */ void addLeg ( const NodeId connectingnode, const int depth, const EdgeLabel edgelabel, LegOccurrences &legoccurrences );
    /* inline */ void addLeftLegs ( Path &path, PathLeg &leg, int &i, Depth olddepth, EdgeLabel lowestlabel, int leftend, int edgesize2 );
    /* inline */ int addLeftLegs ( Path &path, PathLeg &leg, Tuple &tuple, unsigned int legindex, int leftend, int edgesize2 );
    /* inline */ void addRightLegs ( Path &path, PathLeg &leg, int &i, Depth olddepth, EdgeLabel lowestlabel, int rightstart, int nodesize2 );
    /* inline */ int addRightLegs ( Path &path, PathLeg &leg, Tuple &tuple, unsigned int legindex, int rightstart, int nodesize2 );
    PatternTree ( PatternTree &parenttree, unsigned int legindex );
    vector<Tuple> treetuples;
    vector<NodeId> rightmostindexes;
    vector<short> rootpathrelations;
    unsigned int nextprefixindex;
    unsigned int rootpathstart;
    unsigned int nextpathstart;
    unsigned int maxdepth;
    int symmetric; // 0 == not symmetric, 1 == symmetric, even length path, 2 == symmetric, odd length path
    int secondpathleg;
    vector<CloseLegPtr> closelegs;
    friend ostream &operator<< ( ostream &stream, PatternTree &patterntree );
#ifdef GRAPH_OUTPUT
    friend void fillMatrix ( int **A, int &nextnode, int rootnode, NodeLabel rootlabel, 
                  int startpos, int endpos, PatternTree &patterntree );
    NodeLabel tree1rootlabel, tree2rootlabel;
    EdgeLabel rootpathlabel;
#endif
};

#define NONEXTPREFIX ((unsigned int) -1)

#endif
