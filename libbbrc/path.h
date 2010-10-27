// path.h
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

#ifndef PATH_H
#define PATH_H
#include <iostream>
#include <vector>
#include "misc.h"
#include "database.h"
#include "legoccurrence.h"
#include "closeleg.h"
#include "constraints.h"

using namespace std;

struct BbrcPathBbrcTuple {
  BbrcDepth depth;
  BbrcNodeId connectingnode;
  BbrcEdgeLabel edgelabel;
  BbrcNodeLabel nodelabel;
};

struct BbrcPathBbrcLeg {
  BbrcPathBbrcTuple tuple;
  BbrcLegOccurrences occurrences;
};

typedef BbrcPathBbrcLeg *BbrcPathBbrcLegPtr;

class BbrcPath {
  public:
    BbrcPath ( BbrcNodeLabel startnodelabel );
    ~BbrcPath ();
    void expand ();
  private:
    friend class BbrcPatternTree;
    bool is_normal ( BbrcEdgeLabel edgelabel ); // ADDED
    void expand2 (pair<float, string> max);
    BbrcPath ( BbrcPath &parentpath, unsigned int legindex );
    vector<BbrcPathBbrcLegPtr> legs; // pointers used to avoid copy-constructor during a resize of the vector
    vector<BbrcCloseBbrcLegPtr> closelegs;
    vector<BbrcNodeLabel> nodelabels;
    vector<BbrcEdgeLabel> edgelabels;
    int frontsymmetry; // which is lower, the front or front reverse?
    int backsymmetry; // which is lower, the back or back reverse?
    int totalsymmetry; // which is lower, from left to right, or the reverse?

    friend ostream &operator<< ( ostream &stream, BbrcPath &path );
};

#endif
