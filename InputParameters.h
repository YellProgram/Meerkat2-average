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

    InputParameters() :
            reject_outliers(true),
            threshold(3)
    { };
};

InputParameters parse_input(const string& filename);
#endif //MEERKAT2_INPOTPARAMETERS_H
