%module liblast
 
%{
#include "last.h"
%}
 
%include "std_string.i"
using namespace std;
 
%include "std_vector.i"
%template(SVector) std::vector<std::string>;

%include "last.h"
