// misc.h
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

#ifndef MISC_H
#define MISC_H
#include <vector>
#include <stdio.h>
#include <iostream>

using namespace std;

typedef unsigned char EdgeLabel; // combined node-edge label of the input file.
typedef unsigned char NodeLabel;
typedef unsigned short NodeId;
typedef unsigned int Depth; // unsigned int is more efficient than short, but requires more memory...
typedef unsigned int Tid;
typedef unsigned int Frequency;

//extern Frequency minfreq;

#define NOTID ((Tid) -1)
#define NOEDGELABEL ((EdgeLabel) -1)
#define MAXEDGELABEL NOEDGELABEL
#define NONODELABEL ((NodeLabel) -1)
#define NODEPTH ((Depth) -1)
#define NOLEG (-1)
#define NONODE ((NodeId) -1)

// this macro can be used when push_back-ing large structures. In stead of first allocating a local
// variable and pushing this, one first pushes and defines a reference to the space in the vector.
// This avoids re-allocation.
#define vector_push_back(_type,_vector,_var) (_vector).push_back ( _type () ); _type &_var = (_vector).back ();

// can be used to obtain a type when inserting into a map
#define map_insert_pair(_type) typedef typeof(_type) T##_type; pair<T##_type::iterator,bool>

#define store(a,b) { if ( (b).elements.capacity () - (b).elements.size () > (b).elements.size () / 2 ) (a) = (b); else swap ( (a), (b) ); }

#define each(_vector) for (int i = 0 ; i < (int) ( _vector ).size() ; i++ )
#define maxi(a, b) ( (a)>(b) ? (a) : (b) )

extern int level; // 3 : all, 2 : paths and trees, 1 : paths
//extern int maxsize; // maximal size

inline void setmax ( short unsigned int &a, short unsigned int b ) { if ( b > a ) a = b; }

class Statistics {
  public:
    Statistics() : patternsize(0) {}
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



//extern Statistics statistics;


inline vector<string>& operator<<(vector<string>& res, string s) {
  if (s.size()) res.push_back(s);
  return res;
} 



#endif
