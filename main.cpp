#include <iostream>

using namespace std;


int main(int argc, char* argv[]) {
    if(argc < 2 or argc > 2) {
        cout << "usage: meerkat2 filename.mrk" << endl;
        return 0;
    }

    cout << "Meerkat2 v. 0.31" << endl;


    //    ReconstructionParameters par = load_refinement_parameters(argv[1]);


    return 0;
}