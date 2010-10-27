// database.h
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

typedef short InputLastEdgeLabel;
typedef short InputLastNodeLabel;
typedef short InputLastNodeId;
typedef unsigned int LastCombinedInputLabel;

#define combineInputLabels(label1,label2,label3) (label1 | ( ((unsigned int) label2 ) << 16 ) | ( ( (unsigned int) label3 ) << 24 ) )
// maximum 255 node labels for now.

#define NOINPUTEDGELABEL ((InputLastEdgeLabel) -1)
#define NOINPUTNODELABEL ((InputLastNodeLabel) -1)

template<class T>
class Lastpvector {
public:
  T *array;
  int _size;
  Lastpvector<T> ( T *array, int _size ): array ( array ), _size ( _size ) { }
  Lastpvector<T> () { }
  inline int size () const { return _size; }
  void resize ( int s ) { _size = s; }
  void clear () { _size = 0; } // cannot remove allocation, as we are not managing that memory here 
  T &operator[] ( int i ) { return array[i]; }
};

struct LastDatabaseTreeEdge {
  LastEdgeLabel edgelabel;
  LastNodeId tonode;

  LastDatabaseTreeEdge ()  { }

  friend ostream &operator<< ( ostream &stream, LastDatabaseTreeEdge &databasetreeedge );
};

struct LastDatabaseTreeNode {
  LastNodeLabel nodelabel;
  bool incycle;
  Lastpvector<LastDatabaseTreeEdge> edges;

  LastDatabaseTreeNode () { }

  friend ostream &operator<< ( ostream &stream, LastDatabaseTreeNode &databasetreenode );
};

/*
Nodes:			nodes =  [node1, ..., node n]
Edges			edges -> [e1 of n1,...,em of n1, ..., e1 of nn...,ek of nn]
*/

struct LastDatabaseTree {
  LastTid tid, orig_tid;
  int line_nr;
  vector<LastDatabaseTreeNode> nodes;

  LastDatabaseTreeEdge *edges;
  // KS: int activity;
  // KS: float
  float activity;

  // KS: LastDatabaseTree ( LastTid tid , LastTid orig_tid , int line_nr ): tid ( tid ), orig_tid (orig_tid ), line_nr (line_nr), activity ( -1 ) { }
  // KS: initialize to 0.0
  LastDatabaseTree ( LastTid tid , LastTid orig_tid , int line_nr ): tid ( tid ), orig_tid (orig_tid ), line_nr (line_nr), activity ( 0.0 ) { }
  LastDatabaseTree () { }
  
  friend ostream &operator<< ( ostream &stream, LastDatabaseTree &databasetree );
};

typedef LastDatabaseTree *LastDatabaseTreePtr;

struct LastDatabaseLastNodeLabel {
  InputLastNodeLabel inputlabel;
  LastFrequency frequency;
  LastTid lasttid;

  LastLegOccurrences occurrences;
  vector<LastEdgeLabel> frequentedgelabels;

  LastDatabaseLastNodeLabel (): frequency ( 1 ) { }
};

struct LastDatabaseLastEdgeLabel {
  InputLastEdgeLabel inputedgelabel;
  LastNodeLabel tonodelabel, fromnodelabel; 
  LastEdgeLabel edgelabel; // the (order) edge label to which this entry corresponds during the search
  LastFrequency frequency;
  LastTid lasttid;

  LastDatabaseLastEdgeLabel (): frequency ( 1 ) { }
};

class LastDatabase {
  public:
    LastDatabase() {}
    vector<LastDatabaseTreePtr> trees;
    map<LastTid, LastDatabaseTreePtr> trees_map;
    vector<LastDatabaseLastNodeLabel> nodelabels;
    vector<LastDatabaseLastEdgeLabel> edgelabels;
    map<InputLastNodeLabel,LastNodeLabel> nodelabelmap;
    map<LastCombinedInputLabel,LastEdgeLabel> edgelabelmap;
    vector<LastEdgeLabel> edgelabelsindexes; // given an edge label, returns the index of the element in edgelabels in which
    LastEdgeLabel frequentLastEdgeLabelSize () const { return edgelabelsindexes.size (); }
                                         // all information about this edge can be found. Used during the search,
					 // only frequent edge label, node label pairs are stored.

     // NOTE! In the input file, the nodes MUST be listed in pre-order.


     // after "read", determines the frequency of edges, using LastDatabaseLastNodeLabel's edgelasttid/edgelabelfrequency
    void edgecount ();

     // after "edgecount",
     // - removes infrequent data
     // - cleans up the datastructures used until now for counting frequencies
     // - changes the edge label order to optimise the search, fills the database with order numbers instead of
     //   the numbers assigned in the previous Lastlevels; fills edgelabelsindexes.
    void reorder ();

    void printTrees ();
    ~LastDatabase ();
    bool readTreeSmi (string smi, LastTid tid , LastTid orig_tid, int line_nr);
    void readGsp (FILE* input);
    void readTreeGsp (FILE *input, LastTid orig_tid, LastTid tid);
  
  	// Perform DFS through tree to identify cycles
    void determineCycledNodes ( LastDatabaseTreePtr tree, vector<int> &nodestack, vector<bool> &visited1, vector<bool> &visited2 );
};

#endif
