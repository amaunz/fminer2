# Makefile
# © 2008 by Andreas Maunz, andreas@maunz.de, jun 2008

# This file is part of Fminer (fminer).
#
# Fminer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Fminer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Fminer.  If not, see <http://www.gnu.org/licenses/>.

# PATH TO LIBFMINER SOURCE DIRECTORY:

INCLUDE_OB  =  -I/usr/local/include/openbabel-2.0  
INCLUDE_GSL  =   
LDFLAGS_OB  =  -L/usr/local/lib  
LDFLAGS_GSL  =   
INCLUDE_RB  = -I/home/am/.gem/ruby/1.8/gems/rsruby-0.5.1.1/ext

# EXAMPLES:
# ADJUST COMPILER PATH TO OPENBABEL INCLUDE FILES (1st line Linux, 2nd line Windows):
#INCLUDE_OB    = -I/usr/local/include/openbabel-2.0
#INCLUDE_OB  =  -I/usr/local/include/openbabel-2.0  
# ADJUST COMPILER PATH TO GSL INCLUDE FILES (1st line Linux, 2nd line Windows):
#INCLUDE_GSL   = -I/usr/include
#INCLUDE_GSL  =   

# ADJUST LINKER PATH TO OPENBABEL LIBRARY FILES (1st line Linux, 2nd line Windows):
#LDFLAGS_OB    = -L/usr/local/lib
#LDFLAGS_OB  =  -L/usr/local/lib  
# ADJUST LINKER PATH TO GSL LIBRARY FILES (1st line Linux, 2nd line Windows):
#LDFLAGS_GSL   =
#LDFLAGS_GSL  =   

# NORMALLY NO ADJUSTMENT NECESSARY BELOW THIS LINE. Exit and try 'make' now.
# WHAT
PROGRAM = fminer
# OPTIONS
CC            = g++
INCLUDE       = $(INCLUDE_OB) $(INCLUDE_GSL)
LDFLAGS       = $(LDFLAGS_OB) $(LDFLAGS_GSL)
CXXFLAGS      = -Wall -O2 -g
LIBS	      = -lopenbabel -lgsl -lgslcblas

# TARGETS
all: $(PROGRAM) 
$(PROGRAM): main.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) \
	      $(LIBS) \
	      $(LDFLAGS) \
	      -o $@ $<

.PHONY:
clean:
	-rm -rf $(PROGRAM)
