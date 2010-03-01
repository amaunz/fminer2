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



# ADJUST COMPILER PATH TO OPENBABEL INCLUDE FILES (1st line Linux, 2nd line Windows):
INCLUDE_OB  =  -I/usr/local/include/openbabel-2.0  
INCLUDE_OB  =  -I/usr/local/include/openbabel-2.0  
# ADJUST COMPILER PATH TO GSL INCLUDE FILES (1st line Linux, 2nd line Windows):
INCLUDE_GSL  =   
INCLUDE_GSL  =   

# ADJUST LINKER PATH TO OPENBABEL LIBRARY FILES (1st line Linux, 2nd line Windows):
LDFLAGS_OB  =  -L/usr/local/lib  
LDFLAGS_OB  =  -L/usr/local/lib  
# ADJUST LINKER PATH TO GSL LIBRARY FILES (1st line Linux, 2nd line Windows):
LDFLAGS_GSL  =   
LDFLAGS_GSL  =   

# FOR RUBY TARGET: ADJUST COMPILER PATH TO RUBY HEADERS (LINUX)
INCLUDE_RB  = -I/usr/lib/ruby/1.8/i486-linux

# FOR LINUX: INSTALL TARGET DIRECTORY
DESTDIR       = /usr/local/lib/


# NORMALLY NO ADJUSTMENT NECESSARY BELOW THIS LINE. Exit and try 'make' now.
# WHAT
NAME          = last
# OPTIONS
CC            = g++
INCLUDE       = $(INCLUDE_OB) $(INCLUDE_GSL) 
LDFLAGS       = $(LDFLAGS_OB) $(LDFLAGS_GSL)
OBJ           = closeleg.o constraints.o database.o graphstate.o legoccurrence.o path.o patterntree.o last.o
SWIG          = swig
SWIGFLAGS     = -c++ -ruby
ifeq ($(OS), Windows_NT) # assume MinGW/Windows
CXXFLAGS      = -g $(INCLUDE)
LIBS	      = -lm -llibopenbabel-3 -llibgsl -llibgslcblas
LIB1          = lib$(NAME).dll
.PHONY:
all: $(LIB1)
$(LIB1): $(OBJ)
	$(CC) $(LDFLAGS) $(LIBS) -shared -o $@ $^
else                     # assume GNU/Linux
CXXFLAGS      = -g $(INCLUDE) -fPIC 
LIBS_LIB2     = -lopenbabel -lgsl
LIBS          = $(LIBS_LIB2) -ldl -lm -lgslcblas
LIB1          = lib$(NAME).so
LIB1_SONAME   = $(LIB1).1
LIB1_REALNAME = $(LIB1_SONAME).0.1
LIB2          = $(NAME).so
.PHONY:
all: $(LIB1_REALNAME) 
.PHONY:
ruby: $(LIB2)
$(LIB1_REALNAME): $(OBJ)
	$(CC) $(LDFLAGS) $(LIBS) -shared -Wl,-soname,$@ -o $@ $^
	-ln -sf $@ $(LIB1_SONAME)
	-ln -sf $@ $(LIB1)
$(LIB2): $(NAME)_wrap.o $(OBJ)
	$(CC) $(LDFLAGS) -shared $(CXXFLAGS) $^ $(LIBS_LIB2) -o $@
$(NAME)_wrap.o: $(NAME)_wrap.cxx
	$(CC) -c $(CXXFLAGS) $(INCLUDE_RB) $^ -o $@
%.cxx: %.i
	$(SWIG) $(SWIGFLAGS) -o $@ $^
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
	-rm -rf *.o *.cxx $(LIB1) $(LIB1_SONAME) $(LIB1_REALNAME) $(LIB2)
