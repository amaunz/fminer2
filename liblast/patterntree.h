// patterntree.h
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

#ifndef PATTERNTREE_H
#define PATTERNTREE_H
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#include "misc.h"
#include "database.h"
#include "legoccurrence.h"
#include "path.h"
#include "closeleg.h"


using namespace std;

extern bool updated;
extern unsigned int minfreq;

class GSWalk;

struct LastTuple {
  LastDepth depth;
  LastEdgeLabel label;
  LastNodeId connectingnode;
  
  LastTuple ( LastDepth depth, LastEdgeLabel label, LastNodeId connectingnode ): 
    depth ( depth ), label ( label ), connectingnode ( connectingnode ) { }
  LastTuple () { }

  friend bool operator< ( LastTuple &a, LastTuple &b ) { return a.depth > b.depth || ( a.depth == b.depth && a.label < b.label ); }
  friend bool operator== ( LastTuple &a, LastTuple &b ) { return a.depth == b.depth && a.label == b.label; }
  friend ostream &operator<< ( ostream &stream, LastTuple &tuple );
};

struct LastLeg {
  LastTuple tuple;
  LastLastLegOccurrences occurrences;
};

typedef LastLeg *LastLegPtr;

class LastPatternTree {
  public:
    LastPatternTree ( LastPath &path, unsigned int legindex );
    ~LastPatternTree ();
    GSWalk* expand (pair<float, string> max, const int parent_size);
    vector<LastLegPtr> legs; // pointers used to avoid copy-constructor during a resize of the vector
  private:
    void checkIfIndeedNormal ();
    /* inline */ void addExtensionLastLegs ( LastTuple &tuple, LastLastLegOccurrences &legoccurrences );
    /* inline */ void addLastLeg ( const LastNodeId connectingnode, const int depth, const LastEdgeLabel edgelabel, LastLastLegOccurrences &legoccurrences );
    /* inline */ void addLeftLastLegs ( LastPath &path, LastPathLastLeg &leg, int &i, LastDepth olddepth, LastEdgeLabel lowestlabel, int leftend, int edgesize2 );
    /* inline */ int addLeftLastLegs ( LastPath &path, LastPathLastLeg &leg, LastTuple &tuple, unsigned int legindex, int leftend, int edgesize2 );
    /* inline */ void addRightLastLegs ( LastPath &path, LastPathLastLeg &leg, int &i, LastDepth olddepth, LastEdgeLabel lowestlabel, int rightstart, int nodesize2 );
    /* inline */ int addRightLastLegs ( LastPath &path, LastPathLastLeg &leg, LastTuple &tuple, unsigned int legindex, int rightstart, int nodesize2 );
    LastPatternTree ( LastPatternTree &parenttree, unsigned int legindex );
    vector<LastTuple> treetuples;
    vector<LastNodeId> rightmostindexes;
    vector<short> rootpathrelations;
    unsigned int nextprefixindex;
    unsigned int rootpathstart;
    unsigned int nextpathstart;
    unsigned int maxdepth;
    int symmetric; // 0 == not symmetric, 1 == symmetric, even length path, 2 == symmetric, odd length path
    int secondpathleg;
    vector<LastCloseLastLegPtr> closelegs;
    friend ostream &operator<< ( ostream &stream, LastPatternTree &patterntree );
#ifdef GRAPH_OUTPUT
    friend void fillMatrix ( int **A, int &nextnode, int rootnode, LastNodeLabel rootlabel, 
                  int startpos, int endpos, LastPatternTree &patterntree );
    LastNodeLabel tree1rootlabel, tree2rootlabel;
    LastEdgeLabel rootpathlabel;
#endif
};

#define NONEXTPREFIX ((unsigned int) -1)

#endif
