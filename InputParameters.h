//
// Created by Arkadiy on 01/06/2017.
//

#ifndef MEERKAT2_INPOTPARAMETERS_H
#define MEERKAT2_INPOTPARAMETERS_H

#include <string>
#include <vector>

using namespace std;


// the following are a stand alone scripts taken from Meerkat which help to parse a simple input file:
struct ContextAroundPosition {
    string lines_before;
    string current_line;
    string lines_after;
    int pos_within_current_line;
    int line_number;
    string context;
};
vector<string> split(const string& input_string, char where);
ContextAroundPosition get_context(istream& in);
// end of parser

struct InputParameters {
    vector<string> input_files;
    string output_name;
    string symmetry;
    bool reject_outliers;
    double threshold;
    vector<float> slice;
    vector<float> scales;
    bool punch_and_fill;
    float r_punch;
    float r_fill;
    bool bin;
    bool save_reciprocal_space_multipliers;
    bool fft;
    vector<int> binning;
    string centering;
    bool merge_datasets;
    bool trim_for_yell;
    double add_constant;
    bool report_pixel_rint;
    bool report_pixel_variance;
    bool report_pixel_multiplicity;


    InputParameters() :
            reject_outliers(true),
            threshold(3),
            punch_and_fill(false),
            r_punch(0),
            r_fill(0),
            bin(false),
            binning{1,1,1},
            fft(false),
            centering("P"),
            trim_for_yell(false),
            symmetry("1"),
            add_constant(0),
            report_pixel_rint(false),
            report_pixel_variance(false),
            report_pixel_multiplicity(false),
            save_reciprocal_space_multipliers(false),
            merge_datasets(false)
    { };

    bool should_slice() {
        return slice.size()>0;
    }
};

InputParameters parse_input(const string& filename);
#endif //MEERKAT2_INPOTPARAMETERS_H
