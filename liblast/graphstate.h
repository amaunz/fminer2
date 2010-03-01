// graphstate.h
// © 2008 by Andreas Maunz, andreas@maunz.de, jul 2008
// Siegfried Nijssen, snijssen@liacs.nl, feb 2004.

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

#ifndef GRAPHSTATE_H
#define GRAPHSTATE_H
#include <vector>
#include <iostream>
#include <algorithm>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "misc.h"
#include "database.h"
#include "closeleg.h"
#include "patterntree.h"

using namespace std;

typedef unsigned int Mark;

class GSWalk;

class GraphState {
  public:

    struct GSDeletedEdge {
      NodeId tonode, fromnode;
      EdgeLabel edgelabel;
      int postonode, posfromnode;
      bool close;
      Mark cyclemark;
    };

    vector<GSDeletedEdge> deletededges;
    // the current pattern
    vector<Tuple> *treetuples;
    vector<CloseTuple> *closetuples;
    vector<NodeId> nodesinpreorder;

    int backbonelength; // the length of the backbone, in number of nodes
    int startsecondpath; // the position of the second part of the backbone in
                         // the treetuples.
    bool nasty;  // nasty == A-B-A-B-A-B -like cases
    NodeLabel centerlabel;
    EdgeLabel bicenterlabel;
    int closecount;
    bool selfdone; // set by isnormal to store whether the original graph has been
                   // normal-checked; we try to delay this until the very last moment,
                   // as we know that on this graph the normalisation procedure will
                   // have to go through all levels
    
    struct GSEdge {
      int tonode;
      int postonode; // position in the adjacency list of the corresponding reverse edge
      EdgeLabel edgelabel;
      Mark cyclemark;
      bool close; // closing edge
      GSEdge ():cyclemark ( 0 ), close ( false ) { }
      GSEdge ( int tonode, int postonode, EdgeLabel edgelabel, bool close = false ): tonode ( tonode ), postonode ( postonode ), edgelabel ( edgelabel ), cyclemark ( 0 ), close ( close ) { }
    };

    struct GSNode {
      NodeLabel label;
      short unsigned int maxdegree;
      vector<GSEdge> edges;
    };

    //keep for debugging purposes
    void makeState ( DatabaseTree *databasetree );
    void undoState ();
    void insertNode ( NodeLabel nodelabel, short unsigned int maxdegree );
    void deleteNode2 ();
    vector<GSNode> nodes;
    int edgessize;
    short unsigned int getNodeDegree ( int i ) const { return nodes[i].edges.size (); }
    short unsigned int getNodeMaxDegree ( int i ) const { return nodes[i].maxdegree; }
    GraphState  ();
    void determineCycles ( unsigned int usedbit );
    int enumerateSpanning ();
    int is_normal ();
    int normalizetree ();
    int normalizeSelf ();
    void init ();
    void insertStartNode ( NodeLabel nodelabel );
    void deleteStartNode ();
    void insertNode ( int from, EdgeLabel edgelabel, short unsigned int maxdegree );
    void deleteNode ();
    void insertEdge ( int from, int to, EdgeLabel edgelabel );
    void deleteEdge ( int from, int to );
    void deleteEdge ( GSEdge &edge ); // pushes deleted edge on stack
    void reinsertEdge (); // reinserts last edge on the stack
    NodeId lastNode () const { return nodes.size () - 1; }

    void print ( GSWalk* gsw, map<Tid, int> weightmap_a, map<Tid, int> weightmap_i ); 

    void print ( FILE *f );
    void DfsOut(int cur_n, int from_n);
    void to_s ( string& oss );

    void print ( unsigned int frequency );
    void DfsOut(int cur_n, string& oss, int from_n);
    string to_s ( unsigned int frequency );
    string sep();

    void puti(FILE* f, int i);
};

struct GSWNode {
    //      v    <labs>
    // e.g. v    <6 7>
    set<InputNodeLabel> labs;

    int stack(GSWNode n);
    friend ostream& operator<< (ostream &out, GSWNode* n);
};

struct GSWEdge {
    //      e    <to>    <labs>          <occurrences active>    <occurrences inactive>
    // e.g. e    1       <2 1>           {0->2 1->3}             {2 3}
    // meaning                           T. 0 covered by 2 f.
    //                                   on this edge
    int to;
    set<InputEdgeLabel> labs;
    map<Tid, int> a;
    map<Tid, int> i;
    bool deleted;
    int discrete_weight;

    int stack(GSWEdge e);
    static bool lt_to (GSWEdge& e1, GSWEdge& e2){
        if (e1.to < e2.to) return 1;
        return 0;
    }
    static bool equal (GSWEdge* e1, GSWEdge* e2){
        if ((e1->to == e2->to) && std::equal(e1->labs.begin(), e1->labs.end(), e2->labs.begin())) return 1;
        return 0;
    }
    friend ostream& operator<< (ostream &out, GSWEdge* e);
};

class GSWalk {
  public:
      typedef vector<GSWNode> nodevector;      // position represents id, ids are always contiguous
      typedef map<int, map<int, GSWEdge> > edgemap; // key is id of from-node (and to-node)
      nodevector nodewalk;
      edgemap edgewalk;
      nodevector temp_nodewalk;
      edgemap temp_edgewalk;
      vector<int> to_nodes_ex; // nodes that were inserted due to high IDs - must be overwritten
      bool activating;
      int hops;
      float cutoff;
      bool adj_m_sing;
      int adj_m_rank;
      int adj_m_size;

      int conflict_resolution (vector<int> core_ids, GSWalk* s, bool starting=1, int ceiling=0);
        
      static void remove_dups_vector (vector<int>& v) {
          set<int> s;
          pair< set<int>::iterator, bool > result;
          for (vector<int>::iterator i = v.begin(); i < v.end(); ++i) {
              result = s.insert(*i);
              if (!result.second) {
                  v.erase(i);
                  --i; // adjust iterator
              }
          }
      }


      int stack (GSWalk* single, map<int,int> stack_locations);

      void add_edge(int f, GSWEdge e, GSWNode n, bool reorder, vector<int>* core_ids, set<int>* u12);
      void svd();
      void up_edge(int i);
      static bool lt_to_map (pair<int, GSWEdge> a, pair<int, GSWEdge> b) {
        if (a.first < b.first) return 1;
        return 0;
      }
      friend ostream& operator<< (ostream &out, GSWalk* gsw);

      GSWalk() : activating(0), hops(0), cutoff(0.0), adj_m_sing(0), adj_m_rank(0), adj_m_size(0) {
        to_nodes_ex.clear();
      }

};

#endif
