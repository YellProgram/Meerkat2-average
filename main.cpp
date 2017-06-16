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

    cout << "Meerkat-average v. 0.35" << endl;

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
        //TODO: check number of scales is the same as number of input files
        //all necessary parameters are initialised

        for(const auto & s : par.input_files)
            if(!file_exists(s))
                throw FileNotFound(s);

        Slice slice; //By default is everything
        if(par.should_slice())
            slice = Slice(par.slice);

        if(par.scales.size()==0)
            par.scales = vector<float>(par.input_files.size(), 1);
        else if(par.scales.size() != par.input_files.size()){
            cout << "Error: incorrect number of scale coefficients" << endl;
            return 0;
        }

        auto inp = IntensityData<float>::read(par.input_files[0], slice);
        inp.scale(par.scales[0]);
        for(int i = 1; i<par.input_files.size(); ++i) {
            auto t = IntensityData<float>::read(par.input_files[i], slice);
            inp.accumulate(t, par.scales[i]);
        }

        IntensityData<float> res;

        if(par.symmetry == "1")
            res = inp;
        else {
            res = IntensityData<float>::empty(inp);
            average(inp, res, par);
        }

        if(par.punch_and_fill) {
            res.punch_and_fill(par.r_punch, par.r_fill, par.centering);
        }

        if(par.bin)
            res = res.binned(par.binning);

        if(par.fft || par.trim_for_yell)
            res.prepare_for_fft();

        if(par.fft)
            res.FFT();

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
