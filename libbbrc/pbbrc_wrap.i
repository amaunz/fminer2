%module bbrc

%{
#include "bbrc.h"
%}

%include "std_string.i"
using namespace std;

%include "std_vector.i"
%template(SVector) std::vector<std::string>;

%include "bbrc.h"
