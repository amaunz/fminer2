# Makefile
# © 2008 by Andreas Maunz, andreas@maunz.de, jun 2008

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

# ADJUST COMPILER PATH TO OPENBABEL INCLUDE FILES (1st line Linux, 2nd line Windows):
#INCLUDE_OB    = -I/usr/local/include/openbabel-2.0
#INCLUDE_OB   += -I/home/openbabel-2.2.1/include
# ADJUST COMPILER PATH TO GSL INCLUDE FILES (1st line Linux, 2nd line Windows):
#INCLUDE_GSL   = -I/usr/include
#INCLUDE_GSL  += -I/c/Program\ Files/GnuWin32/include/

# ADJUST LINKER PATH TO OPENBABEL LIBRARY FILES (1st line Linux, 2nd line Windows):
#LDFLAGS_OB    = -L/usr/local/lib
#LDFLAGS_OB   += -L/c/Users/am/Desktop/fminer
# ADJUST LINKER PATH TO GSL LIBRARY FILES (1st line Linux, 2nd line Windows):
#LDFLAGS_GSL   =
#LDFLAGS_GSL  += -L/c/Program\ Files/GnuWin32/bin/

# FOR RUBY TARGET: ADJUST COMPILER PATH TO RUBY HEADERS (LINUX)
#INCLUDE_RB    = -I/usr/lib/ruby/1.8/i486-linux/
#INCLUDE_JAVA  = -I/usr/lib/jvm/java-6-sun-1.6.0.10/include/ -I/usr/lib/jvm/java-6-sun-1.6.0.10/include/linux/

# FOR LINUX: INSTALL TARGET DIRECTORY
DESTDIR       = /usr/local/lib/


# NORMALLY NO ADJUSTMENT NECESSARY BELOW THIS LINE. Exit and try 'make' now.
# WHAT
NAME          = bbrc
# OPTIONS
CC            = g++
INCLUDE       = $(INCLUDE_OB) $(INCLUDE_GSL) 
LDFLAGS       = $(LDFLAGS_OB) $(LDFLAGS_GSL)
OBJ           = closeleg.o constraints.o database.o graphstate.o legoccurrence.o path.o patterntree.o bbrc.o
SWIG          = swig
SWIGFLAGS     = -c++
ifeq ($(OS), Windows_NT) # assume MinGW/Windows
CXXFLAGS      = -O3 -g $(INCLUDE)
LIBS	      = -lm -llibopenbabel-3 -llibgsl -llibgslcblas
LIB1          = lib$(NAME).dll
.PHONY:
all: $(LIB1)
$(LIB1): $(OBJ)
	$(CC) $(LDFLAGS) $(LIBS) -shared -o $@ $^
else                     # assume GNU/Linux
CXXFLAGS      = -O3 -g $(INCLUDE) -fPIC
LIBS_LIB2     = -lopenbabel -lgsl
LIBS          = $(LIBS_LIB2) -ldl -lm -lgslcblas
LIB1          = lib$(NAME).so
LIB1_SONAME   = $(LIB1).1
LIB1_REALNAME = $(LIB1_SONAME).0.1
LIB2          = $(NAME).so
LIB3          = lib$(NAME).so
.PHONY:
all: $(LIB1_REALNAME) 
.PHONY:
ruby: $(LIB2)
java: $(LIB3)
$(LIB1_REALNAME): $(OBJ)
	$(CC) $(LDFLAGS) $(LIBS) -shared -Wl,-soname,$@ -o $@ $^
	-ln -sf $@ $(LIB1_SONAME)
	-ln -sf $@ $(LIB1)

$(LIB2): clean $(OBJ) r$(NAME)_wrap.i
	$(SWIG) $(SWIGFLAGS) -ruby -o $(NAME)_wrap.cxx r$(NAME)_wrap.i
	$(CC) -c $(CXXFLAGS) $(INCLUDE_RB) $(NAME)_wrap.cxx -o $(NAME)_wrap.o
	$(CC) $(LDFLAGS) -shared $(CXXFLAGS) $(OBJ) $(NAME)_wrap.o $(LIBS_LIB2) -o $@

$(LIB3): clean $(OBJ) j$(NAME)_wrap.i
	$(SWIG) $(SWIGFLAGS) -java -o $(NAME)_wrap.cxx j$(NAME)_wrap.i
	$(CC) -c $(CXXFLAGS) $(INCLUDE_JAVA) $(NAME)_wrap.cxx -o $(NAME)_wrap.o
	$(CC) $(LDFLAGS) -shared $(CXXFLAGS) $(OBJ) $(NAME)_wrap.o $(LIBS_LIB2) -o $@

endif

install: $(LIB1_REALNAME)
	cp -P $(LIB1)* $(DESTDIR)

# FILE TARGETS
.o: .cpp.h
	$(CC) -Wall $(CXXFLAGS) $(LIBS) $@
# HELPER TARGETS
.PHONY:
doc: Doxyfile Mainpage.h *.h
	-doxygen $<
.PHONY:
clean:
	-rm -rf *.o *.cxx libfminer.java  Fminer.*  libfminerJNI.*  SVector.*  SWIGTYPE_p_ChisqConstraint.java  SWIGTYPE_p_FILE.* $(LIB1) $(LIB1_SONAME) $(LIB1_REALNAME) $(LIB2)
