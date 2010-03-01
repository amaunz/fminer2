// database.cpp
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

#include "database.h"
#include "constraints.h"
#include <algorithm>
#include <iostream>

namespace fm {
    extern bool aromatic;
    extern unsigned int minfreq;
}

ostream &operator<< ( ostream &stream, DatabaseTreeEdge &databasetreeedge ) {
  stream << "DatabaseTreeEdge; edgelabel: " << databasetreeedge.edgelabel << "; tonode: " << databasetreeedge.tonode << endl;
  return stream;
}

ostream &operator<< ( ostream &stream, DatabaseTreeNode &databasetreenode ) {
  stream << "DatabaseTreeNode; label: " << databasetreenode.nodelabel << "; edges: " << endl;
  for (int i = 0; i < databasetreenode.edges.size (); i++ )
    stream << databasetreenode.edges[i];
  stream << endl;
  return stream;
}

ostream &operator<< ( ostream &stream, DatabaseTree &databasetree ) {
  stream << "DatabaseTree; tid: " << databasetree.tid << "; nodes: " << endl;
  for (unsigned int i = 0; i < databasetree.nodes.size (); i++ )
    stream << databasetree.nodes[i];
  stream << endl;
  return stream;
}





bool Database::readTreeSmi (string smi, Tid tid, Tid orig_tid, int line_nr) {

    OBMol mol;

    istringstream iss (smi, istringstream::in);
    OBConversion conv(&iss,&cout);

    // read the molecule
    conv.SetInAndOutFormats("SMI","SDF");

    if (!conv.ReadString(&mol,smi)) {
        cerr << "Error during conversion" << endl;
        return(0);
    }   

    if (!mol.DeleteHydrogens()) {
        cerr << "Unable to delete hydrogens" << endl;
        return(0);
    }

    // create and store new tree object
    DatabaseTreePtr tree = new DatabaseTree ( tid , orig_tid , line_nr );


    // SEGFAULT
    trees.push_back(tree);
    trees_map[orig_tid] = tree;

    int nodessize = 0, edgessize = 0;
    static vector<DatabaseTreeNode> nodes;
    static vector<vector<DatabaseTreeEdge> > edges;
    nodes.resize ( 0 );

//    cerr << "Atoms are (Type(ID)):" << endl;
    OBAtomIterator atom;
    mol.BeginAtom(atom);

	///////////
	// NODES //
	///////////

//    cerr << endl;
    do {

        //cerr << " " << (*atom)->GetType() << " (idx " << (*atom)->GetIdx() << ")" << endl;

        InputNodeLabel inputnodelabel=0;

        // set atom type as label
        // code for 'c' is set to -1 (aromatic carbon).
        if (fm::aromatic) {
            ((*atom)->IsAromatic() && ((*atom)->GetAtomicNum()==6)) ? inputnodelabel = 254 : inputnodelabel = (*atom)->GetAtomicNum();
        }
        else inputnodelabel = (*atom)->GetAtomicNum();
        nodessize++;

        // Insert into map, using subsequent numbering for internal labels:
    	// node nr. 1, node nr. 2, ...
	    // Store direction inputlabel -> internal label
    	//cerr << "NodeLabelMap: insert " << inputnodelabel << " --> " << nodelabels.size() << endl;
        map_insert_pair ( nodelabelmap ) p = nodelabelmap.insert ( make_pair ( inputnodelabel, (nodelabels.size()) ) );

        // Store direction internal label -> node
        // if node has NOT been present, label it and set frequency to 1
        if ( p.second ) {
          vector_push_back ( DatabaseNodeLabel, nodelabels, nodelabel );
          nodelabel.inputlabel = inputnodelabel;
          nodelabel.occurrences.parent = NULL;
          nodelabel.occurrences.number = 1;
          nodelabel.lasttid = tid;
          //cerr << "Created node label " << nodelabel.inputlabel << " (freq " << nodelabel.frequency <<  ")" << endl;
        }
        // if node has been present, increase frequency
        else {
          DatabaseNodeLabel &nodelabel = nodelabels[p.first->second];
          if ( nodelabel.lasttid != tid ) nodelabel.frequency++;
          nodelabel.lasttid = tid;
          //cerr << "Updated node label " << nodelabel.inputlabel << " (freq " << nodelabel.frequency << ")" << endl;
        }

        // Tree node
        vector_push_back ( DatabaseTreeNode, nodes, node );
        node.nodelabel = p.first->second;                           // refer to nodelabel
    //    node.atom = (*atom);                                        // attach OB atom
        node.incycle = false;
        //cerr << "Created tree node for OB index " << node.atom->GetIdx() << " (nodelabel " << (int) node.nodelabel << ", nodes[] size " << nodessize << ")" << endl;

    } while (mol.NextAtom(atom));


    // copy nodes to tree and prepare edges storage size
    // edges[nodeid] gives the edges going out of node with id 'nodeid'
    tree->nodes.reserve ( nodessize );
    if ( edges.size () < (unsigned int) nodessize )
        edges.resize ( nodessize );					// edges stored per node
    for ( int i = 0; i < nodessize; i++ ) {
        edges[i].resize ( 0 );						// no edges yet
        tree->nodes.push_back ( nodes[i] );
    }
    //cerr << endl;

    InputEdgeLabel inputedgelabel;
    InputNodeId nodeid1, nodeid2;

    //cerr << "Bonds are (idx[label] bo idx[label]):" << endl;
    OBBondIterator bond;

    ///////////
    // EDGES //
    ///////////
    
    if (mol.BeginBond(bond)) {

//        cerr << endl;
        
        do {

            nodeid1 = (InputNodeId) ((*bond)->GetBeginAtomIdx())-1;     // USE OB INDICES (same as nodelabel+1)!
            nodeid2 = (InputNodeId) ((*bond)->GetEndAtomIdx())-1;       //
            int bondorder = (*bond)->GetBondOrder();

            // set input edge label
            inputedgelabel = bondorder;
            if (fm::aromatic && (*bond)->IsAromatic()) inputedgelabel = 4;

//            cerr << nodeid1 << inputedgelabel << "(" << (*bond)->IsAromatic() << ")" << nodeid2 << " ";
            NodeLabel node1label = tree->nodes[nodeid1].nodelabel;
            NodeLabel node2label = tree->nodes[nodeid2].nodelabel;
            
            //cerr << "(" << (*bond)->GetBeginAtomIdx() << "[" << (int) node1label << "] " << (*bond)->GetBondOrder() << " " << (*bond)->GetEndAtomIdx() << "[" << (int) node2label << "]"<< ")" << endl;

            
            // Direction of edge always from low to high
            if ( node1label > node2label ) {
                NodeLabel temp = node1label;
                node1label = node2label;
                node2label = temp;
            }

            // Create combined input label for edge
            CombinedInputLabel combinedinputlabel;
            combinedinputlabel = combineInputLabels ( inputedgelabel, node1label, node2label );

            // Insert into map, analoguous to nodes, using subsequent numbering for internal labels:
            // edge nr. 1, edge nr. 2, ...
            
            // Direction inputlabel -> internal label
            //cerr << "EdgeLabelMap: insert " << combinedinputlabel << "-->" << edgelabels.size() << endl;
            map_insert_pair ( edgelabelmap ) p = edgelabelmap.insert ( make_pair ( combinedinputlabel, edgelabels.size () ) );
     
            // Direction internal label -> edge
            if ( p.second ) {
              vector_push_back ( DatabaseEdgeLabel, edgelabels, edgelabel );
              edgelabel.fromnodelabel = node1label;	// directed edges
              edgelabel.tonodelabel = node2label;
              edgelabel.inputedgelabel = inputedgelabel;
              edgelabel.lasttid = tid;
              //cerr << "Created edge " << edgelabel.inputedgelabel << " (" << (int) edgelabel.fromnodelabel << "-->" << (int) edgelabel.tonodelabel << ")" << ":" << edgelabel.frequency << endl;
            }
            else {
              DatabaseEdgeLabel &edgelabel = edgelabels[p.first->second];
              if ( edgelabel.lasttid != tid )
                edgelabel.frequency++;
                edgelabel.lasttid = tid;
                //cerr << "Updated edge " << edgelabel.inputedgelabel << " (" << (int) edgelabel.fromnodelabel << "-->" << (int) edgelabel.tonodelabel << ")" << ":" << edgelabel.frequency << endl;
            }
        
            // Tree edge (2 versions)
            vector_push_back ( DatabaseTreeEdge, edges[nodeid1], edge );
            edge.edgelabel = p.first->second;
            edge.tonode = nodeid2;
    //        edge.bond = (*bond);
                
            vector_push_back ( DatabaseTreeEdge, edges[nodeid2], edge2 );
            edge2.edgelabel = p.first->second;
            edge2.tonode = nodeid1;
    //        edge2.bond = (*bond);
           
            edgessize++;

        } while (mol.NextBond(bond));
    }

    // copy edges to tree
    tree->edges = new DatabaseTreeEdge[edgessize * 2];	// allocate space
                                                        // two instances per edge (both directions)
    int pos = 0;
    for ( int i = 0; i < nodessize; i++ ) {				// for all nodes...
        int s = edges[i].size ();
        tree->nodes[i].edges._size = s;						// edges stored in s steps from ...
        tree->nodes[i].edges.array = tree->edges + pos;		// ... offset for that node's edges
        for ( int j = 0; j < s; j++, pos++ ) {
            //cerr << i << " " << j << " L: " << edges[i][j].tonode << endl;
            tree->edges[pos] = edges[i][j];					// flattened storage
        }
    }
  
    ////////////
    // CYCLES //
    ////////////

    static vector<int> nodestack;
    static vector<bool> visited1, visited2;
    nodestack.resize ( 0 );
    visited1.resize ( 0 );
    visited1.resize ( nodessize, false );
    visited2.resize ( 0 );
    visited2.resize ( nodessize, false );
  


    // for every node...
    for ( int i = 0; i < nodessize; i++ ) {
        if ( !visited1[i] ) {
            nodestack.push_back ( i );
            visited1[i] = visited2[i] = true;		// visited1: has been reached by DFS
											        // visited2: is on the current path (indicator for cycle)
	        // ...perform DFS
            determineCycledNodes ( tree, nodestack, visited1, visited2 );
            visited2[i] = false;
            nodestack.pop_back ();
        }
    }
    return(1);
}

void Database::readGsp (FILE* input) {
  Tid tid2 = 0; 

  char array[100];
  fgets ( array, 100, input );
  string tree_s = array; Tid orig_tid = (unsigned int) atoi((tree_s.substr(tree_s.find_first_of("123456789"))).c_str());

  while ( !feof ( input ) ) {
    readTreeGsp ( input, tid2, orig_tid );
    fgets ( array, 100, input );
    tree_s = array; orig_tid = (unsigned int) atoi((tree_s.substr(tree_s.find_first_of("123456789"))).c_str());
    tid2++;
  }

}

int readint ( FILE *input ) {
  char car = fgetc ( input );
  while ( car < '0' || car > '9' ) {
    if ( feof ( input ) )
      return -1;
    car = fgetc ( input );
  }
  int n = car - '0';
  car = fgetc ( input );
  while ( car >= '0' && car <= '9' ) {
    n = n * 10 + car - '0';
    car = fgetc ( input );
  }

  return n;
}


char readcommand ( FILE *input ) {
  char car = fgetc ( input );
  while ( car < 'a' || car > 'z' ) {
    if ( feof ( input ) ) {
      return -1;
    }
    car = fgetc ( input );
  }
  return car;
}

void Database::readTreeGsp ( FILE *input, Tid tid , Tid orig_tid) {
  InputNodeLabel inputnodelabel;

  DatabaseTreePtr tree = new DatabaseTree ( tid , orig_tid , tid );
  trees.push_back ( tree );
  trees_map[orig_tid] = tree;

  char command;
  int dummy;
  int nodessize = 0, edgessize = 0;
  command = readcommand ( input );
  
  static vector<DatabaseTreeNode> nodes;
  static vector<vector<DatabaseTreeEdge> > edges;
  nodes.resize ( 0 );

  while ( command == 'v' ) {
    dummy = readint ( input );
    inputnodelabel = readint ( input );
    if ( dummy != nodessize ) {
      cerr << "Error reading input file - node number does not correspond to its position." << endl;
      exit ( 1 );
    }
    nodessize++;

    map_insert_pair ( nodelabelmap ) p = nodelabelmap.insert ( make_pair ( inputnodelabel, nodelabels.size () ) );
    if ( p.second ) {
      vector_push_back ( DatabaseNodeLabel, nodelabels, nodelabel );
      nodelabel.inputlabel = inputnodelabel;
      nodelabel.occurrences.parent = NULL;
      nodelabel.occurrences.number = 1;
      nodelabel.lasttid = tid;
    }
    else {
      DatabaseNodeLabel &nodelabel = nodelabels[p.first->second];
      if ( nodelabel.lasttid != tid )
        nodelabel.frequency++;
      nodelabel.lasttid = tid;
    }
    
    vector_push_back ( DatabaseTreeNode, nodes, node );
    node.nodelabel = p.first->second;
    node.incycle = false;

    command = readcommand ( input );
  }

  tree->nodes.reserve ( nodessize );
  if ( (int) edges.size () < nodessize )
    edges.resize ( nodessize );
  for ( int i = 0; i < nodessize; i++ ) {
    edges[i].resize ( 0 );
    tree->nodes.push_back ( nodes[i] );
  }
  
  InputEdgeLabel inputedgelabel;
  InputNodeId nodeid1, nodeid2;

  while ( !feof ( input ) && command == 'e' ) {
    nodeid1 = readint ( input );
    nodeid2 = readint ( input );
    inputedgelabel = readint ( input );
    NodeLabel node2label = tree->nodes[nodeid2].nodelabel;
    NodeLabel node1label = tree->nodes[nodeid1].nodelabel;
    CombinedInputLabel combinedinputlabel;
    if ( node1label > node2label ) {
      NodeLabel temp = node1label;
      node1label = node2label;
      node2label = temp;
    }
    combinedinputlabel = combineInputLabels ( inputedgelabel, node1label, node2label );

    map_insert_pair ( edgelabelmap ) p = edgelabelmap.insert ( make_pair ( combinedinputlabel, edgelabels.size () ) );
    if ( p.second ) {
      vector_push_back ( DatabaseEdgeLabel, edgelabels, edgelabel );
      edgelabel.fromnodelabel = node1label;
      edgelabel.tonodelabel = node2label;
      edgelabel.inputedgelabel = inputedgelabel;
      edgelabel.lasttid = tid;
    }
    else {
      DatabaseEdgeLabel &edgelabel = edgelabels[p.first->second];
      if ( edgelabel.lasttid != tid )
        edgelabel.frequency++;
      edgelabel.lasttid = tid;
    }

    vector_push_back ( DatabaseTreeEdge, edges[nodeid1], edge );
    edge.edgelabel = p.first->second;
    edge.tonode = nodeid2;

    vector_push_back ( DatabaseTreeEdge, edges[nodeid2], edge2 );
    edge2.edgelabel = p.first->second;
    edge2.tonode = nodeid1;

    edgessize++;

    command = readcommand ( input );
    if (command == 't') fseek (input, -1, SEEK_CUR);
  }

  tree->edges = new DatabaseTreeEdge[edgessize * 2];
  int pos = 0;
  for ( int i = 0; i < nodessize; i++ ) {
    int s = edges[i].size ();
    tree->nodes[i].edges._size = s;
    tree->nodes[i].edges.array = tree->edges + pos;
    for ( int j = 0; j < s; j++, pos++ ) {
      tree->edges[pos] = edges[i][j];
    }
  }
  
  static vector<int> nodestack;
  static vector<bool> visited1, visited2;
  nodestack.resize ( 0 );
  visited1.resize ( 0 );
  visited1.resize ( nodessize, false );
  visited2.resize ( 0 );
  visited2.resize ( nodessize, false );
  for ( int i = 0; i < nodessize; i++ ) {
    if ( !visited1[i] ) {
      nodestack.push_back ( i );
      visited1[i] = visited2[i] = true;
      determineCycledNodes ( tree, nodestack, visited1, visited2 );
      visited2[i] = false;
      nodestack.pop_back ();
    }
  }
}









// TREE                      
// determine internal cycles 

void Database::determineCycledNodes ( DatabaseTreePtr tree, vector<int> &nodestack, vector<bool> &visited1, vector<bool> &visited2 ) {
    int node = nodestack.back ();
    pvector<DatabaseTreeEdge> &edges = tree->nodes[node].edges;	// jump to beginning of (flat) edge storage

    for ( int i = 0; i < edges.size (); i++ ) {		            // edges.size() answered by pvector
        if ( !visited1[edges[i].tonode] ) {				        // [] operator = array + i steps!!
            //cerr << edges[i].tonode << " (" << nodestack.size() << ")" << endl;
            nodestack.push_back ( edges[i].tonode );
          
            visited1[edges[i].tonode] = visited2[edges[i].tonode] = true;
            determineCycledNodes ( tree, nodestack, visited1, visited2 );
            nodestack.pop_back ();
            visited2[edges[i].tonode] = false;
        }
        else {                                  // exclude other direction of the same edge................................//
            if ( visited2[edges[i].tonode] && ( nodestack.size () == 1 || nodestack[nodestack.size () - 2] != edges[i].tonode ) ) {
		        //cerr << "Detected cycle: " << edges[i].tonode << endl;
                int j = nodestack.size () - 1;
                while ( nodestack[j] != edges[i].tonode ) {
		            //cerr << "Marking: " << nodestack[j] << " (" << (tree->nodes[nodestack[j]]).atom->GetType() << ")" << endl;
                    tree->nodes[nodestack[j]].incycle = true;
                    j--;
                }
		        //cerr << "Marking: " << nodestack[j] << " (" << (tree->nodes[nodestack[j]]).atom->GetType() << ")" << endl;
                tree->nodes[nodestack[j]].incycle = true;
            }
        }
    }
}

void Database::edgecount () {
  for (unsigned int i = 0; i < edgelabels.size (); i++ ) {                              // DATABASE                    
    if ( edgelabels[i].frequency >= fm::minfreq ) {                                         // if edge is frequent...      
      nodelabels[edgelabels[i].tonodelabel].frequentedgelabels.push_back ( i );         // ... store it at the to-node 
      if ( edgelabels[i].fromnodelabel != edgelabels[i].tonodelabel )                   // ... and also (if different) 
        nodelabels[edgelabels[i].fromnodelabel].frequentedgelabels.push_back ( i );     // ... at the from-node        
    }
  }
}

class EdgeLabelsIndexesSort:public std::binary_function<int,int,bool> {
    const vector<DatabaseEdgeLabel> &edgelabels;
  public:
    EdgeLabelsIndexesSort ( const vector<DatabaseEdgeLabel> &edgelabels ) : edgelabels ( edgelabels ) { }
    bool operator () ( int a, int b ) const {
      return edgelabels[a].frequency < edgelabels[b].frequency;
    }
};

bool operator< ( const DatabaseTreeEdge &a, const DatabaseTreeEdge &b ) {
  return a.edgelabel < b.edgelabel;
}

void Database::reorder () {
    

    // PHASE I: LABEL EDGES ACCORDING TO FREQUENCY
  //  cerr << "REORDER" << endl;


    // gather frequent edgelabels and sort according to frequency
    edgelabelsindexes.reserve ( edgelabels.size () );
    for (unsigned int i = 0; i < edgelabels.size (); i++ ) {
        if ( edgelabels[i].frequency >= fm::minfreq )
            edgelabelsindexes.push_back ( i );                                                              
    }

  //each(edgelabelsindexes) {
  //      cerr << (int) edgelabelsindexes[i] << " (" 
  //         << nodelabels[edgelabels[edgelabelsindexes[i]].fromnodelabel].inputlabel
  //         << " => "
  //         << nodelabels[edgelabels[edgelabelsindexes[i]].tonodelabel].inputlabel
  //         << ")" << endl;
  //}
  //cerr << endl;

    sort ( edgelabelsindexes.begin (), edgelabelsindexes.end (), EdgeLabelsIndexesSort ( edgelabels ) );

  //each(edgelabelsindexes) {
  //    cerr << (int) edgelabelsindexes[i] << " (" 
  //         << nodelabels[edgelabels[edgelabelsindexes[i]].fromnodelabel].inputlabel
  //         << " => "
  //         << nodelabels[edgelabels[edgelabelsindexes[i]].tonodelabel].inputlabel
  //         << ")" << endl;
  //}
  //cerr << endl;


    // Now, use the ranked indices to re-number frequent edge labels with their rank
    for (unsigned int i = 0; i < edgelabelsindexes.size (); i++ ) {
        edgelabels[edgelabelsindexes[i]].edgelabel = i;                 // fill in the edge labels for the first time
//     #define DEBUG
     #ifdef DEBUG
    DatabaseEdgeLabel &label = edgelabels[edgelabelsindexes[i]];
    cout << (int) nodelabels[label.tonodelabel].inputlabel 
         << "[" << (int) label.inputedgelabel << "]" 
	 << (int) nodelabels[label.fromnodelabel].inputlabel <<"-->" << i <<endl;
     #endif                                                                   // the edgelabel is the rank found by REORDER
    }


    // PHASE II: REMOVE INFREQUENT EDGES FROM NODES
//    cerr << "REMOVE" << endl;

    for ( Tid i = 0; i < trees.size (); i++ ) {
        DatabaseTree &tree = * (trees[i]);                                          // for every tree i...
        for ( NodeId j = 0; j < tree.nodes.size (); j++ ) {                         // for every node j...
  //        cerr << endl;
            DatabaseTreeNode &node = tree.nodes[j];
            if ( nodelabels[node.nodelabel].frequency >= fm::minfreq ) {                  // ...check its frequency...
                DatabaseNodeLabel &nodelabel = nodelabels[node.nodelabel];

  //            cerr << "Leg Occurence for node " << nodelabel.inputlabel
  //                 << ": " << tree.tid << " " << nodelabel.occurrences.elements.size() << " " << j << " " << NONODE << endl;
                nodelabel.occurrences.elements.push_back ( LegOccurrence ( tree.tid, (OccurrenceId) nodelabel.occurrences.elements.size (), j, NONODE )  );
                                                                                        // ...and push occurence in database
                int k = 0;
  //            cerr << "node " << (int) node.nodelabel  << " (" << nodelabels[node.nodelabel].inputlabel << ")"  << endl;
                for ( int l = 0; l < node.edges.size (); l++ ) {                        // For each edge l going out of j...
                    
                                        
                    EdgeLabel lab = node.edges[l].edgelabel;                            // ... (with label lab)...
                    if ( edgelabels[lab].frequency >= fm::minfreq ) {                       // ... check its frequency...

  //                    DatabaseTreeEdge& edge = node.edges[l];
  //                    cerr << "  edge " << (int) edge.edgelabel << " moved from " << l << " to " << k << endl;

                        node.edges[k].edgelabel = edgelabels[lab].edgelabel;            // ... and overwrite old edges
                        node.edges[k].tonode = node.edges[l].tonode;
                        //node.edges[k].bond = node.edges[l].bond;
                        k++;
                    }
                }
  //            cerr << "Truncating " << node.edges.size() - k << " edges" << endl;
                node.edges.resize ( k );                                                // truncate the rest
            }
            else node.edges.clear ();                                                   // truncate all edges
        }
    }


}

void Database::printTrees () {
  for (unsigned int i = 0; i < trees.size (); i++ )
    cout << trees[i];
}

Database::~Database () {
  for (unsigned int i = 0; i < trees.size (); i++ )
    delete trees[i];
}
