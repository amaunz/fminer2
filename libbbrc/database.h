// database.h
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

#ifndef DATABASE_H
#define DATABASE_H
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdio>

#include <openbabel/mol.h>
#include "openbabel/obconversion.h"

#include "legoccurrence.h"
#include "misc.h"

using namespace std;
using namespace OpenBabel;

typedef short InputBbrcEdgeLabel;
typedef short InputBbrcNodeLabel;
typedef short InputBbrcNodeId;
typedef unsigned int BbrcCombinedInputLabel;

#define combineInputLabels(label1,label2,label3) (label1 | ( ((unsigned int) label2 ) << 16 ) | ( ( (unsigned int) label3 ) << 24 ) )
// maximum 255 node labels for now.

#define NOINPUTEDGELABEL ((InputBbrcEdgeLabel) -1)
#define NOINPUTNODELABEL ((InputBbrcNodeLabel) -1)

template<class T>
class Bbrcpvector {
public:
  T *array;
  int _size;
  Bbrcpvector<T> ( T *array, int _size ): array ( array ), _size ( _size ) { }
  Bbrcpvector<T> () { }
  inline int size () const { return _size; }
  void resize ( int s ) { _size = s; }
  void clear () { _size = 0; } // cannot remove allocation, as we are not managing that memory here 
  T &operator[] ( int i ) { return array[i]; }
};

struct BbrcDatabaseTreeEdge {
  BbrcEdgeLabel edgelabel;
  BbrcNodeId tonode;

  BbrcDatabaseTreeEdge ()  { }

  friend ostream &operator<< ( ostream &stream, BbrcDatabaseTreeEdge &databasetreeedge );
};

struct BbrcDatabaseTreeNode {
  BbrcNodeLabel nodelabel;
  bool incycle;
  Bbrcpvector<BbrcDatabaseTreeEdge> edges;

  BbrcDatabaseTreeNode () { }

  friend ostream &operator<< ( ostream &stream, BbrcDatabaseTreeNode &databasetreenode );
};

/*
Nodes:			nodes =  [node1, ..., node n]
Edges			edges -> [e1 of n1,...,em of n1, ..., e1 of nn...,ek of nn]
*/

struct BbrcDatabaseTree {
  BbrcTid tid, orig_tid;
  int line_nr;
  vector<BbrcDatabaseTreeNode> nodes;

  BbrcDatabaseTreeEdge *edges;
  // KS: int activity;
  // KS: float
  float activity;
  float weight;

  // KS: BbrcDatabaseTree ( BbrcTid tid , BbrcTid orig_tid , int line_nr ): tid ( tid ), orig_tid (orig_tid ), line_nr (line_nr), activity ( -1 ) { }
  // KS: initialize to 0.0
  BbrcDatabaseTree ( BbrcTid tid , BbrcTid orig_tid , int line_nr ): tid ( tid ), orig_tid (orig_tid ), line_nr (line_nr), activity ( 0.0 ), weight (0.0) { }
  BbrcDatabaseTree () { }
  
  friend ostream &operator<< ( ostream &stream, BbrcDatabaseTree &databasetree );
};

typedef BbrcDatabaseTree *BbrcDatabaseTreePtr;

struct BbrcDatabaseBbrcNodeLabel {
  InputBbrcNodeLabel inputlabel;
  BbrcFrequency frequency;
  BbrcTid lasttid;

  BbrcLegOccurrences occurrences;
  vector<BbrcEdgeLabel> frequentedgelabels;

  BbrcDatabaseBbrcNodeLabel (): frequency ( 1 ) { }
};

struct BbrcDatabaseBbrcEdgeLabel {
  InputBbrcEdgeLabel inputedgelabel;
  BbrcNodeLabel tonodelabel, fromnodelabel; 
  BbrcEdgeLabel edgelabel; // the (order) edge label to which this entry corresponds during the search
  BbrcFrequency frequency;
  BbrcTid lasttid;

  BbrcDatabaseBbrcEdgeLabel (): frequency ( 1 ) { }
};

class BbrcDatabase {
  public:
    BbrcDatabase() {}
    vector<BbrcDatabaseTreePtr> trees;
    map<BbrcTid, BbrcDatabaseTreePtr> trees_map;
    vector<BbrcDatabaseBbrcNodeLabel> nodelabels;
    vector<BbrcDatabaseBbrcEdgeLabel> edgelabels;
    map<InputBbrcNodeLabel,BbrcNodeLabel> nodelabelmap;
    map<BbrcCombinedInputLabel,BbrcEdgeLabel> edgelabelmap;
    vector<BbrcEdgeLabel> edgelabelsindexes; // given an edge label, returns the index of the element in edgelabels in which
    BbrcEdgeLabel frequentBbrcEdgeLabelSize () const { return edgelabelsindexes.size (); }
                                         // all information about this edge can be found. Used during the search,
					 // only frequent edge label, node label pairs are stored.

     // NOTE! In the input file, the nodes MUST be listed in pre-order.


     // after "read", determines the frequency of edges, using BbrcDatabaseBbrcNodeLabel's edgelasttid/edgelabelfrequency
    void edgecount ();

     // after "edgecount",
     // - removes infrequent data
     // - cleans up the datastructures used until now for counting frequencies
     // - changes the edge label order to optimise the search, fills the database with order numbers instead of
     //   the numbers assigned in the previous Bbrclevels; fills edgelabelsindexes.
    void reorder ();

    void printTrees ();
    ~BbrcDatabase ();
    bool readTreeSmi (string smi, BbrcTid tid , BbrcTid orig_tid, int line_nr);
    void readGsp (FILE* input);
    void readTreeGsp (FILE *input, BbrcTid orig_tid, BbrcTid tid);
  
  	// Perform DFS through tree to identify cycles
    void determineCycledNodes ( BbrcDatabaseTreePtr tree, vector<int> &nodestack, vector<bool> &visited1, vector<bool> &visited2 );
};

#endif
