// path.h
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

#ifndef PATH_H
#define PATH_H
#include <iostream>
#include <vector>
#include <assert.h>
#include "misc.h"
#include "database.h"
#include "legoccurrence.h"
#include "closeleg.h"
#include "constraints.h"

using namespace std;

class GSWalk;

struct LastPathLastTuple {
  LastDepth depth;
  LastNodeId connectingnode;
  LastEdgeLabel edgelabel;
  LastNodeLabel nodelabel;
};

struct LastPathLastLeg {
  LastPathLastTuple tuple;
  LastLegOccurrences occurrences;
};

typedef LastPathLastLeg *LastPathLastLegPtr;

class LastPath {
  public:
    LastPath ( LastNodeLabel startnodelabel );
    ~LastPath ();
    void expand ();
  private:
    friend class LastPatternTree;
    bool is_normal ( LastEdgeLabel edgelabel ); // ADDED
    GSWalk* expand2 (pair<float, string> max, const int parent_size);
    LastPath ( LastPath &parentpath, unsigned int legindex );
    vector<LastPathLastLegPtr> legs; // pointers used to avoid copy-constructor during a resize of the vector
    vector<LastCloseLastLegPtr> closelegs;
    vector<LastNodeLabel> nodelabels;
    vector<LastEdgeLabel> edgelabels;
    int frontsymmetry; // which is lower, the front or front reverse?
    int backsymmetry; // which is lower, the back or back reverse?
    int totalsymmetry; // which is lower, from left to right, or the reverse?

    friend ostream &operator<< ( ostream &stream, LastPath &path );
};

#endif
