%module libfminer

%{
#include "fminer.h"
%}

%include "std_string.i"
typedef std::string String;
using namespace std;

%include "std_vector.i"
%template(SVector) std::vector<std::string>;

%include "fminer.h"
