// path.h
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

struct PathTuple {
  Depth depth;
  NodeId connectingnode;
  EdgeLabel edgelabel;
  NodeLabel nodelabel;
};

struct PathLeg {
  PathTuple tuple;
  LegOccurrences occurrences;
};

typedef PathLeg *PathLegPtr;

class Path {
  public:
    Path ( NodeLabel startnodelabel );
    ~Path ();
    void expand ();
  private:
    friend class PatternTree;
    bool is_normal ( EdgeLabel edgelabel ); // ADDED
    void expand2 (pair<float, string> max);
    Path ( Path &parentpath, unsigned int legindex );
    vector<PathLegPtr> legs; // pointers used to avoid copy-constructor during a resize of the vector
    vector<CloseLegPtr> closelegs;
    vector<NodeLabel> nodelabels;
    vector<EdgeLabel> edgelabels;
    int frontsymmetry; // which is lower, the front or front reverse?
    int backsymmetry; // which is lower, the back or back reverse?
    int totalsymmetry; // which is lower, from left to right, or the reverse?

    friend ostream &operator<< ( ostream &stream, Path &path );
};

#endif
