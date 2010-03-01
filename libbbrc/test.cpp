#include "fminer.h"
#include <iostream>
#include <string.h>
using namespace std;

Fminer* MyFminer;
int main(int argc, char *argv[], char *envp[]) {
    MyFminer= new Fminer();
    MyFminer->AddCompound ("COC1=CC=C(C=C1)C2=NC(=C([NH]2)C3=CC=CC=C3)C4=CC=CC=C4" , 1);
    MyFminer->AddCompound ("O=C1NC(=S)NC(=O)C1C(=O)NC2=CC=CC=C2" , 2);
    // ... continue adding compounds
    MyFminer->AddActivity((bool) true, 1);
    MyFminer->AddActivity((bool) false, 2);
    // ... continue adding activities (true for active, false for inactive)
    cerr << MyFminer->GetNoCompounds() << " compounds" << endl;
    // Toy example: special settings for mining all fragments
    MyFminer->SetChisqSig(0); // use no significance constraint
    MyFminer->SetRefineSingles(true); // refine structures with support 1
    // gather results for every root node in vector instead of immediate output
    MyFminer->SetConsoleOut(false);
    for ( int j = 0; j < (int) MyFminer->GetNoRootNodes(); j++ ) {
        vector<string>* result = MyFminer->MineRoot(j);
        for ( int i = 0; i < result->size(); i++) {
            cout << (*result)[i] << endl;
        }
    }
    delete MyFminer;
}

