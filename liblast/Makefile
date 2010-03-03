# Makefile
# (c) 2010 by Andreas Maunz, andreas@maunz.de, feb 2010

# This file is part of LibFminer (libfminer).
#
# LibFminer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# LibFminer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with LibFminer.  If not, see <http://www.gnu.org/licenses/>.

INCLUDE_OB  =  -I/usr/local/include/openbabel-2.0  
INCLUDE_GSL  =   
LDFLAGS_OB  =  -L/usr/local/lib  
LDFLAGS_GSL  =   
INCLUDE_RB  = -I/usr/lib/ruby/1.8/i486-linux
RUBY = ruby1.8

# ADJUST COMPILER PATH TO OPENBABEL INCLUDE FILES:
#INCLUDE_OB    = -I/usr/local/include/openbabel-2.0
# ADJUST COMPILER PATH TO GSL INCLUDE FILES:
#INCLUDE_GSL   = -I/usr/include

# ADJUST LINKER PATH TO OPENBABEL LIBRARY FILES:
#LDFLAGS_OB    = -L/usr/local/lib
# ADJUST LINKER PATH TO GSL LIBRARY FILES:
#LDFLAGS_GSL   =

# FOR RUBY TARGET: ADJUST COMPILER PATH TO RUBY HEADERS (LINUX)
#INCLUDE_RB    = -I/usr/lib/ruby/1.8/i486-linux/
#RUBY          = ruby1.8

# NORMALLY NO ADJUSTMENT NECESSARY BELOW THIS LINE. Exit and try 'make' now.
# WHAT
NAME          = last
# OPTIONS
CC            = g++
INCLUDE       = $(INCLUDE_OB) $(INCLUDE_GSL) 
LDFLAGS       = $(LDFLAGS_OB) $(LDFLAGS_GSL)
OBJ           = closeleg.o constraints.o database.o graphstate.o legoccurrence.o path.o patterntree.o last.o
CXXFLAGS      = -g $(INCLUDE) -fPIC 
LIBS_LIB2     = -lopenbabel -lgsl
LIBS          = $(LIBS_LIB2) -ldl -lm -lgslcblas
LIB1          = lib$(NAME).so
LIB1_SONAME   = $(LIB1).2
LIB1_REALNAME = $(LIB1_SONAME).0
LIB2          = $(NAME).so

# FILES
$(LIB1_REALNAME): $(OBJ)
	$(CC) $(LDFLAGS) $(LIBS) -shared -Wl,-soname,$@ -o $@ $^
	-ln -sf $@ $(LIB1_SONAME)
	-ln -sf $@ $(LIB1)

$(LIB2): $(OBJ) $(NAME)_wrap.cxx
	$(CC) -c $(CXXFLAGS) $(INCLUDE_RB) $(NAME)_wrap.cxx -o $(NAME)_wrap.o
	$(CC) $(LDFLAGS) -shared $(CXXFLAGS) $(OBJ) $(NAME)_wrap.o $(LIBS) -o $@

.o: .cpp.h
	$(CC) -Wall $(CXXFLAGS) $(LIBS) $@

ctest: test.cpp
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(OBJ) $(LIBS) $< -o $@


# TESTS + HELPERS
.PHONY:
all: $(LIB1_REALNAME) 
.PHONY:
ruby: $(LIB2)
.PHONY:
cpptest: ctest
	./$<
.PHONY:
rbtest: test.rb
	$(RUBY) $<
.PHONY:
doc: Doxyfile Mainpage.h *.h
	-doxygen $<
.PHONY:
clean:
	-rm -rf *.o $(LIB1) $(LIB1_SONAME) $(LIB1_REALNAME) $(LIB2) ctest
