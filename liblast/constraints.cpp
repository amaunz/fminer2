// constraints.cpp
// © 2008 by Andreas Maunz, andreas@maunz.de, jun 2008

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

#include "constraints.h"

float ChisqConstraint::ChiSq(float x, float y, bool decide_activating) {

        float ea = 0.0, ei = 0.0, impact = 0.0;
        
        impact = x/(float)n;
        ea = na * impact; 
        ei = ni * impact; 

        if (decide_activating) {
            if (y>ea) activating=1; else activating=0;
        }

        if (ea>0 && ei>0) chisq = (y-ea-0.5)*(y-ea-0.5)/ea + (x-y-ei-0.5)*(x-y-ei-0.5)/ei;

        return(chisq);

}
