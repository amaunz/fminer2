#include "last.h"
#include <iostream>
#include <string.h>
using namespace std;

Last* MyFminer;
int main(int argc, char *argv[], char *envp[]) {
    MyFminer= new Last();
    // Adjust settings
    MyFminer->SetMaxHops(25);
    MyFminer->SetConsoleOut(false);
    // Add compounds below. IMPORTANT! DO NOT CHANGE SETTINGS AFTER ADDING COMPOUNDS!
    MyFminer->AddCompound ("O=C(C(C(C=C3)=CC=C3O)=CO2)C1=C2C=C(O)C=C1O" , 1);
    MyFminer->AddCompound ("Oc1ccc(cc1)[C@@H]2Cc3ccc(O)cc3OC2" , 2);
    MyFminer->AddCompound ("O=C1C(C3=CC=C(O)C=C3)=COC2=C1C=CC(O)=C2" , 3);
    MyFminer->AddCompound ("O=C1C(C3=CC=C(OC)C=C3)=COC2=C1C=CC(O)=C2" , 4);
    MyFminer->AddCompound ("OC1=CC=C(CCCCCCCC)C=C1" , 5);
    MyFminer->AddCompound ("C1(C=CC=CC=1C(=C(Cl)Cl)C2=CC=C(C=C2)Cl)Cl" , 6);
    MyFminer->AddCompound ("O=C(C1=C(C=CC=C1)C(=O)OCC(CCCC)CC)OCC(CCCC)CC" , 7);
    MyFminer->AddCompound ("Oc1cc(O)cc2CCCCC[C@@H](O)CCC[C@H](C)OC(=O)c12" , 8);
    MyFminer->AddCompound ("O=C1C2=C(C=C(C=C2O)O)OC(=C1O)C3=CC(=C(C=C3)O)O" , 9);
    MyFminer->AddCompound ("C1(=C(C(=O)C2=C(O1)C=C(C=C2)O)O)C3=CC(O)=C(C=C3)O" , 10);
    // ... continue adding compounds
    MyFminer->AddActivity((bool) 1.0, 1);
    MyFminer->AddActivity((bool) 1.0, 2);
    MyFminer->AddActivity((bool) 1.0, 3);
    MyFminer->AddActivity((bool) 1.0, 4);
    MyFminer->AddActivity((bool) 1.0, 5);
    MyFminer->AddActivity((bool) 1.0, 6);
    MyFminer->AddActivity((bool) 1.0, 7);
    MyFminer->AddActivity((bool) 1.0, 8);
    MyFminer->AddActivity((bool) 0.0, 9);
    MyFminer->AddActivity((bool) 0.0, 10);
    // ... continue adding activities (true for active, false for inactive)
    cerr << MyFminer->GetNoCompounds() << " compounds" << endl;
    // gather results for every root node in vector instead of immediate output
    for ( int j = 0; j < (int) MyFminer->GetNoRootNodes(); j++ ) {
        vector<string>* result = MyFminer->MineRoot(j);
        for ( int i = 0; i < result->size(); i++) {
            cout << (*result)[i] << endl;
        }
    }
    delete MyFminer;
}

