// fminer.h
// © 2010 by Andreas Maunz, andreas@maunz.de, feb 2010

/*
    This file is part of Fminer (fminer). 

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

#ifndef FMINER_H
#define FMINER_H

#include <vector>
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <errno.h>
#include <stdio.h>
#include <iomanip>
#include <algorithm>

#define each(_vector) for (int i = 0 ; i < (int) ( _vector ).size() ; i++ )

typedef unsigned char EdgeLabel; // combined node-edge label of the input file.
typedef unsigned char NodeLabel;
typedef unsigned short NodeId;
typedef unsigned int Depth; // unsigned int is more efficient than short, but requires more memory...
typedef unsigned int Tid;

class Fminer {

public:

    /** @name Inits
     *  Initializer functions.
     */
    //@{
    Fminer () {}
    Fminer (int _type, float _minfreq) {}
    Fminer (int _type, float _minfreq, float chisq_val, bool _do_backbone) {}
    virtual ~Fminer() {}
    //@}
    /** @name Getters
     *  Getter functions.
     */
    //@{
    virtual bool GetConsoleOut() = 0;
    virtual bool GetRegression() = 0;
    //@}
    /** @name Setters
     *  Setter functions.
     */
    //@{
    virtual void SetMinfreq(float val) = 0;
    virtual bool SetType(int val) = 0;
    virtual bool SetBackbone(bool val) = 0;
    virtual bool SetDynamicUpperBound(bool val) = 0;
    virtual bool SetPruning(bool val) = 0;
    virtual bool SetConsoleOut(bool val) = 0;
    virtual void SetAromatic(bool val) = 0;
    virtual bool SetRefineSingles(bool val) = 0;
    virtual void SetDoOutput(bool val) = 0;
    virtual bool SetBbrcSep(bool val) = 0;
    virtual bool SetChisqActive(bool val) = 0;
    virtual bool SetChisqSig(float _chisq_val) = 0;
    virtual bool SetRegression(bool val) = 0;
    virtual bool SetMaxHops(int val) = 0;
    //@}
    /** @name Others
     *  Other functions.
     */
    //@{
    virtual std::vector<std::string>* MineRoot(unsigned int j) = 0;
    virtual void ReadGsp(FILE* gsp) = 0;
    virtual bool AddCompound(std::string smiles, unsigned int comp_id) = 0;
    virtual bool AddActivity(float act, unsigned int comp_id) = 0;
    virtual int GetNoRootNodes() = 0;
    virtual int GetNoCompounds() = 0;
    //@}

};

// the types of the class factories
typedef Fminer* create0_t();
typedef Fminer* create2_t(int _type, float _minfreq);
typedef Fminer* create4_t(int _type, float _minfreq, float _chisq_val, bool _do_backbone);
typedef void destroy_t(Fminer*);
typedef void usage_f();

#endif
