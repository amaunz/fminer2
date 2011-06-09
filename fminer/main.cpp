// main.cpp
// © 2010 by Andreas Maunz, andreas@maunz.de, jun 2010

/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <getopt.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <cassert>
#include <dlfcn.h>

#include "fminer.h"

using namespace std;

class Statistics;
extern Statistics* statistics;

Fminer* fminer;
destroy_t* destroy_lib;
void* Lib;


// helper routines
void puti ( FILE *f, int i ) {
  char array[100]; 
  int k = 0;
  do {
    array[k] = ( i % 10 ) + '0';
    i /= 10;
    k++; 
  }
  while ( i != 0 ); 
  do {
    k--;
    putc ( array[k], f );
  } while ( k );
}

void remove_dos_cr(string* str) {
	string nl = "\r"; 
	for (string::size_type i = str->find(nl); i!=string::npos; i=str->find(nl)) str->erase(i,1); // erase dos cr
}

void read_gsp (char* graph_file) {
    FILE *input = fopen (graph_file, "r");
    if (!input) {
        cerr << "Error opening file '" << graph_file << "': " << strerror(errno) << "." << endl;
        exit(1);
    }
    fminer->ReadGsp(input);
}

void read_smi (char* graph_file) {
    Tid tree_id = 0;

    ifstream input;
    string line;
    int line_nr=0;
    string tmp_field="";

    // Open input stream
    input.open(graph_file);
    if (!input) {
        cerr << "Error opening file '" << graph_file << "': " << strerror(errno) << "." << endl;
        exit(1);
    }

    // draw line from input
    int field_nr = 2; // initialize to 'valid' pattern
    while (getline(input, line)) {
        line_nr++;
        istringstream iss(line); // line as iss
        string id = "";

        // if previous line has had two entries, i.e. format was correct, then re-initialize
        if (field_nr==2) field_nr=0;
        else {
            cerr << "Error! Line no. " << line_nr-1 << " in file '" << graph_file << "' is not in correct input format. Please refer to README." << endl;
            exit(1);
        }

        while(getline(iss, tmp_field, '\t')) {  // split at tabs
            if (field_nr == 0) {        // ID
                tree_id = (unsigned int) atoi (tmp_field.c_str());
                if (tree_id == 0) { 
                    cerr << "Error! Invalid ID: '" << tmp_field << "' in file  '" << graph_file << "', line " << line_nr+1 << "." << endl; 
                    exit(1); 
                }
                field_nr=1;
            }
            else if (field_nr == 1) {   // SMILES
                fminer->AddCompound (tmp_field , tree_id);
                field_nr=2;
            }
			else {
                cerr << "Error! Line no. " << line_nr << " in file '" << graph_file << "' is not in correct input format. Please refer to README." << endl;
				exit(1);
			}
        }
    }
    input.close();

}

void read_act (char* act_file, bool regr) {
    ifstream input;
    string line;
    string tmp_field;
    string act_name;
    Tid tid=0;
    unsigned int line_nr = 0;
    
    // Open input stream
    input.open(act_file);
    if (!input) {
        cerr << "Error opening file '" << act_file << "': " << strerror(errno) << "." << endl;
        exit(1);
    }

    // draw line from input
    unsigned int field_nr = 3; // initialize to 'valid' pattern
	while (getline(input, line)) {
        line_nr++;
		istringstream iss(line); // line as iss

        if (field_nr==3) field_nr=0;
        else {
            cerr << "Error! Line no. " << line_nr-1 << " in file '" << act_file << "' is not in correct input format. Please refer to README." << endl;
            exit(1);
        }

		while(getline(iss, tmp_field, '\t')) {	// split at tabs
			remove_dos_cr(&tmp_field);

			if (field_nr == 0) {		    // ID
                tid = (Tid) atoi(tmp_field.c_str());
                if (tid == 0) { 
                    cerr << "Error! Invalid ID: '" << tmp_field << "' in file '" << act_file << "', line " << line_nr << "." << endl; 
                    exit(1); 
                }
                field_nr=1;
			}

			else if (field_nr == 1) {	    // ACTIVITY NAME
				act_name = tmp_field;
                field_nr=2;
			}

			else if (field_nr == 2) {	    // ACTIVITY VALUES
                stringstream str;

                // KS: str  << tmp_field; int act_value; str >> act_value;
                // KS: convert to float
                str  << tmp_field; float act_value; str >> act_value;
                
                /* KS:
                if ((act_value != 0) && act_value != 1) { 
                    cerr << "Error! Invalid activity: '" << tmp_field << "' in file '" << act_file << "', line " << line_nr << "." << endl; 
                    exit(1);
                }
                */
                // KS: check float for equality

                // KS: fminer->AddActivity((bool) act_value, tid);
                // KS: Do not convert to bool
                fminer->AddActivity(act_value, tid);
                
                field_nr=3;
			}

			else {
                cerr << "Error! Line no. " << line_nr << " in file '" << act_file << "' is not in correct input format. Please refer to README." << endl;
				exit(1);
			}

		}

        if (field_nr != 3) {
            cerr << "Error! Line no. " << line_nr << " in file '" << act_file << "' is not in correct input format. Please refer to README." << endl;
            exit(1);
        }


	}


}


// main
int main(int argc, char *argv[], char *envp[]) {


    Frequency def_minfreq = 2;
    Frequency minfreq = def_minfreq;
    bool arg_minfreq = 0;

    int def_type = 2;
    int type = def_type;
    bool arg_type = 0; // arg_* needed to also capture argument with default val

    float def_chisq_sig = 0.95;
    float chisq_sig = def_chisq_sig;
    bool arg_chisq_sig = 0;

    bool def_refine_singles = false;
    bool refine_singles = def_refine_singles;

    bool def_aromatic = true;
    bool aromatic = def_aromatic;

    bool def_do_pruning = true;
    bool do_pruning = def_do_pruning;

    bool def_do_backbone = true;
    bool do_backbone = def_do_backbone;

    bool def_adjust_ub = true;
    bool adjust_ub = def_adjust_ub;

    bool def_do_output = true;
    bool do_output = def_do_output;

    bool def_do_regression = false;
    bool do_regression = def_do_regression;
 
    int def_max_hops = 1000;
    int max_hops = def_max_hops;
    bool arg_max_hops = 0;
   
    int status=1;
    const char* program_name = argv[0];
    char* graph_file = NULL;
    char* act_file = NULL;
    char* lib_path = NULL;

    
    // FILE ARGUMENT READ: STATUS 1
	if (argc>1) {
	    if (argv[1][0]!='-') {
            lib_path = argv[1]; //set lib path
            if (argc>3) {
               if (argv[argc-2][0]!='-') {
                   graph_file = argv[argc-2]; 
                   status=0;
                   if (argv[argc-1][0]=='-') {
                       status=1;
                   }
                   else {
                       act_file = argv[argc-1];
                   }
               }
               else {
                   if (argv[argc-1][0]=='-') {
                       status=1;
                   }
                   else {
                      graph_file = argv[argc-1];
                      status = 0;
                   }
               }
            }
            else if (argc==3){
               if (argv[argc-1][0]=='-') {
                   status=1;
               }
               else {
                   graph_file = argv[argc-1];
                   status = 0;
               }
            }
            else status=1;
        }
        else status=1;
	}
	else status=1;


    // OPTIONS ARGUMENT READ
    char c;
    const char* const short_options = "f:l:p:saubdogm:h";
    const struct option long_options[] = {
        {"minfreq",                1, NULL, 'f'},
        {"level",                  1, NULL, 'l'},
        {"p-value",                1, NULL, 'p'},
        {"refine-singles",         0, NULL, 's'},
        {"no-aromaticity",         0, NULL, 'a'},
        {"no-upper-bound-pruning", 0, NULL, 'u'},
        {"no-bbr-classes",         0, NULL, 'b'},
        {"no-dynamic-ub",          0, NULL, 'd'},
        {"no-output",              0, NULL, 'o'},
        {"regression",             0, NULL, 'g'},
        {"max-hops",               1, NULL, 'm'},
        {"help",                   0, NULL, 'h'},
        {NULL,                     0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
        switch(c) {
        case 'f':
            minfreq = atoi(optarg);
            arg_minfreq = 1;
            break;
        case 'l':
            type = atoi (optarg);
            arg_type = 1;
            break;
        case 'p':
            chisq_sig = atof (optarg);
            arg_chisq_sig = 1;
            if (!act_file) status = 1;
            break;
        case 's':
            refine_singles = true;
            break;
        case 'a':
            aromatic = false;
            break;
        case 'u':
            do_pruning = false;
            if (!act_file) status = 1;
            break;
        case 'b':
            do_backbone = false;
            if (!act_file) status = 1;
            break;
        case 'd':
            adjust_ub = false;
            if (!act_file) status = 1;
            break;
        case 'o':
            do_output = false;
            break;
        case 'g':
            do_regression = true;
            if (!act_file) status = 1;
            break;
        case 'm':
            max_hops = atoi(optarg);
            arg_max_hops = 1;
            break;
        case 'h':
            if ((argc>1) && (argv[1][0]!='-')) status=2;
            break;
        case '?':
            status=1;
            if ((argc>1) && (argv[1][0]!='-')) status=2;
            break;
        default: 
            abort();
        }
    }


    if (status == 0) {
        // INTEGRITY CONSTRAINTS AND HELP OUTPUT
        //  ----------- !du ---------      ----------- !db ---------
        if ((adjust_ub && !do_pruning) || (adjust_ub && !do_backbone)) status = 2; 
        if (do_regression && (!adjust_ub || !do_backbone || !do_pruning) ) status = 2; // KS: enforce d,b,u flags not set
    }

    bool input_smi = false, input_gsp = false;
    string graph_file_str;
    if (graph_file) {
        graph_file_str = graph_file;
        string graph_file_suffix = graph_file_str.substr(graph_file_str.find_last_of("."));
        if (graph_file_suffix == ".smi") { input_smi=true; }
        else if (graph_file_suffix == ".gsp") { input_gsp=true; }
        else { cerr << "Suffix " << graph_file_suffix << " unknown!" << endl; status=1;}
    }

    if (status==0 || status == 2) {

       // Check which library
       Lib = dlopen(lib_path, RTLD_LAZY);
       if (!Lib) {
            cerr << "Cannot load library: " << dlerror() << '\n';
            return 1;
       }


       else { 
            create0_t* create_lib = (create0_t*) dlsym(Lib, "create0");
            const char* dlsym_error = dlerror();
            if (dlsym_error) {
                cerr << "Cannot load symbol create: " << dlsym_error << '\n';
                return 1;
            }
            destroy_lib = (destroy_t*) dlsym(Lib, "destroy");
            dlsym_error = dlerror();
            if (dlsym_error) {
                cerr << "Cannot load symbol destroy: " << dlsym_error << '\n';
                return 1;
            }    

            bool all_args_good = 1; // switched to 0 in case of illegal arguments, as reported by lib.

            if (graph_file && act_file) {
                fminer = create_lib();
                if (do_regression != def_do_regression) all_args_good *= fminer->SetRegression(do_regression);
                if (type != def_type || arg_type) all_args_good *= fminer->SetType(type);
                if (minfreq != def_minfreq || arg_minfreq) fminer->SetMinfreq(minfreq);
                if (chisq_sig != def_chisq_sig || arg_chisq_sig) all_args_good *= fminer->SetChisqSig(chisq_sig);
                if (refine_singles != def_refine_singles) all_args_good *= fminer->SetRefineSingles(refine_singles);
                if (aromatic != def_aromatic) fminer->SetAromatic(aromatic);
                if (adjust_ub != def_adjust_ub) all_args_good *= fminer->SetDynamicUpperBound(adjust_ub);
                if (do_pruning != def_do_pruning) all_args_good *= fminer->SetPruning(do_pruning);
                if (do_backbone != def_do_backbone) all_args_good *= fminer->SetBackbone(do_backbone);
                if (do_output != def_do_output) fminer->SetDoOutput(do_output);
                //if (bbrc_sep != def_bbrc_sep) all_args_good *= fminer->SetBbrcSep(bbrc_sep); // Disabled for console output. Set manually to true and disable console output.
                if (max_hops != def_max_hops || arg_max_hops)  all_args_good *= fminer->SetMaxHops(max_hops);
            }

            else if (graph_file) {
                fminer = create_lib();
                if (type != def_type || arg_type) all_args_good *= fminer->SetType(type);
                if (minfreq != def_minfreq || arg_minfreq) fminer->SetMinfreq(minfreq);
                if (refine_singles != def_refine_singles) all_args_good *= fminer->SetRefineSingles(refine_singles);
                if (aromatic != def_aromatic) fminer->SetAromatic(aromatic);
                if (do_output != def_do_output) fminer->SetDoOutput(do_output);
                fminer->SetChisqActive(false);
            }

            if (!all_args_good) status = 2;
        }
    }

    if (status > 0) {
        cerr << endl;
        cerr << "Fminer v2.0, Andreas Maunz, 2010" << endl;
        cerr << "Usage 1: " << program_name << " <Library> <Options> <Graphs> <Activities>" << endl;
        cerr << "Usage 2: " << program_name << " <Library> <Options> <Graphs>" << endl;
        cerr << endl;
        cerr << "File formats:" << endl;
        cerr << "       <Library>    Plug-in library to use (/path/to/libbbrc.so or /path/to/liblast.so)." << endl;
        cerr << "       <Graphs>     File should have suffix .smi or .gsp, indicating SMILES or gSpan format." << endl;
        cerr << "       <Activities> File must be in Activity format (suffix not relevant)." << endl;
        cerr << endl;

    }

    if (status==1) {
        cerr << "Use '" << program_name << " <Library> -h' for additional information." << endl;
        cerr << endl;
    }

    else if (status==2) {
        usage_f* print_usage = (usage_f*) dlsym(Lib, "usage");
        print_usage();
        cerr << "See README for additional information." << endl;
        cerr << endl;
        destroy_lib(fminer); 
        dlclose(Lib); 
        return 1;
    }

    if (status > 0) return 1;
   

 
    // status 0 -> go ahead
    fminer->SetConsoleOut(true);
    
    //////////
    // READ //
    //////////


    if (graph_file && act_file) {
        cerr << "Reading compounds..." << endl;
        if (input_smi) read_smi (graph_file);
        else if (input_gsp) read_gsp (graph_file);
        cerr << "Reading activities..." << endl;
        read_act (act_file, fminer->GetRegression());
    }
    
    else if (graph_file) {
        cerr << "Reading compounds..." << endl;
        if (input_smi) read_smi(graph_file);
        else if (input_gsp) read_gsp(graph_file);
    }

    //////////
    // MINE //
    //////////
    
    cerr << "Mining fragments..." << endl;

    cerr << fminer->GetNoCompounds() << " compounds" << endl;
    clock_t t1 = clock ();
    for ( int j = 0; j < (int) fminer->GetNoRootNodes(); j++ ) {
        vector<string>* result = fminer->MineRoot(j);
        if (!fminer->GetConsoleOut()) { 
            each (*result) {
                cout << (*result)[i] << endl;
            }
        }
    }
    clock_t t2 = clock ();
//  statistics->print();
    cerr << "Approximate total runtime: " << ( (float) t2 - t1 ) / CLOCKS_PER_SEC << "s" << endl;
    
    destroy_lib(fminer);
    dlclose(Lib);

}
