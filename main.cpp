#include <iostream>
#include <fstream>
#include "InputParameters.h"
#include "Exceptions.h"
#include "IntensityData.h"
#include "Symmetry.h"

using namespace std;


bool file_exists(const string& filename) {
    ifstream in(filename);
    return in.good();
}

int main(int argc, char* argv[]) {
    if(argc < 2 or argc > 2) {
        cout << "usage: meerkat-average filename.mrk" << endl;
        return 0;
    }

    cout << "Meerkat-average v. 0.31" << endl;

    //    ReconstructionParameters par = load_refinement_parameters(argv[1]);

    try {
        InputParameters par = parse_input(argv[1]);

        cout << "Averaging datasets: " << endl;
        for(const auto& s : par.input_files)
            cout << "    " << s << endl;

        cout << "Output name: " << par.output_name << endl;

        if(file_exists(par.output_name)) {
            cout << "Warning: file \"" << par.output_name << "\" will be overwritten" << endl;
        }

        cout << endl;

        //CHECKS HERE or in parser:
        //all nesessary parameters are initialised

        for(const auto & s : par.input_files)
            if(!file_exists(s))
                throw FileNotFound(s);

        //MAIN LOGIC HERE
        auto inp = IntensityData<float>::read(par.input_files[0]);
        auto res = IntensityData<float>::empty(inp);
        average(inp, res, par);
        res.save(par.output_name);



    }catch(const FileNotFound& f_err) {
        cout << "Error: file \"" << f_err.filename << "\" does not exist." << endl;
        return 0;
    } catch (const std::bad_alloc&) {
        cout  << endl << "Error: not enough operating memory." << endl;
        return 0;
    } catch (const ParserError& e) {
        cout << e.description;
        return 0;
    }

    return 0;
}
