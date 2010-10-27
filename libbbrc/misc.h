// misc.h
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

#ifndef MISC_H
#define MISC_H
#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>

using namespace std;

typedef unsigned char BbrcEdgeLabel; // combined node-edge label of the input file.
typedef unsigned char BbrcNodeLabel;
typedef unsigned short BbrcNodeId;
typedef unsigned int BbrcDepth; // unsigned int is more efficient than short, but requires more memory...
typedef unsigned int BbrcTid;
typedef unsigned int BbrcFrequency;

//extern BbrcFrequency minfreq;

#define NOTID ((BbrcTid) -1)
#define NOEDGELABEL ((BbrcEdgeLabel) -1)
#define MAXEDGELABEL NOEDGELABEL
#define NONODELABEL ((BbrcNodeLabel) -1)
#define NODEPTH ((BbrcDepth) -1)
#define NOLEG (-1)
#define NONODE ((BbrcNodeId) -1)

// this macro can be used when push_back-ing large structures. In stead of first allocating a local
// variable and pushing this, one first pushes and defines a reference to the space in the vector.
// This avoids re-allocation.
#define vector_push_back(_type,_vector,_var) (_vector).push_back ( _type () ); _type &_var = (_vector).back ();

// can be used to obtain a type when inserting into a map
#define map_insert_pair(_type) typedef typeof(_type) T##_type; pair<T##_type::iterator,bool>

#define store(a,b) { if ( (b).elements.capacity () - (b).elements.size () > (b).elements.size () / 2 ) (a) = (b); else swap ( (a), (b) ); }

#define each(_vector) for (int i = 0 ; i < (int) ( _vector ).size() ; i++ )
#define maxi(a, b) ( (a)>(b) ? (a) : (b) )

extern int Bbrclevel; // 3 : all, 2 : paths and trees, 1 : paths
//extern int maxsize; // maximal size

inline void Bbrcsetmax ( short unsigned int &a, short unsigned int b ) { if ( b > a ) a = b; }

class BbrcStatistics {
  public:
    BbrcStatistics() : patternsize(0) {}
    vector<unsigned int> frequenttreenumbers;
    vector<unsigned int> frequentpathnumbers;
    vector<unsigned int> frequentgraphnumbers;
    int patternsize;
    void print () {
        int total = 0, total2 = 0, total3 = 0;
        for (unsigned int i = 0; i < frequenttreenumbers.size (); i++ ) {
          cerr << "Frequent " << i + 2
               << " cyclic graphs: " << frequentgraphnumbers[i]
               << " real trees: " << frequenttreenumbers[i]
               << " paths: " << frequentpathnumbers[i]
               << " total: " << frequentgraphnumbers[i] + frequenttreenumbers[i] + frequentpathnumbers[i] << endl;
          total += frequentgraphnumbers[i];
          total2 += frequenttreenumbers[i];
          total3 += frequentpathnumbers[i];
        }
        cerr << "TOTAL:" << endl
           << "Frequent cyclic graphs: " << total << " real trees: " << total2 << " paths: " << total3 << " total: " << total + total2 + total3 << endl;
    }  
};



//extern BbrcStatistics statistics;


inline vector<string>& operator<<(vector<string>& res, string s) {
  if (s.size()) res.push_back(s);
  return res;
} 



#endif
