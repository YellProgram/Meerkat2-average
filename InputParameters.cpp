//
// Created by Arkadiy on 01/06/2017.
//

#include "InputParameters.h"
#include "Exceptions.h"
#include <string>
#include <fstream>
#include <sstream>
#include <set>

using namespace std;

template<class T>
bool isIn(const T& target, const std::set<T>& the_set)
{
    return the_set.find(target) != the_set.end();
}

bool operator== (const ContextAroundPosition& a, const ContextAroundPosition& b) {
    return a.lines_after==b.lines_after && a.lines_before==b.lines_before && a.current_line==b.current_line && a.pos_within_current_line==b.pos_within_current_line;
}

ContextAroundPosition get_context(istream& in) {
    in.clear();

    if(in.eof())
    {
        in.clear();
        in.seekg(0,ios_base::end);
    }

    int context_position = in.tellg();
    vector<string> lines;

    in.seekg(0,ios_base::beg);

    string line;
    if(context_position==0) { //in case the line is empty
        getline(in,line);
        lines.push_back(line);
    } else
        while(int(in.tellg()) < context_position and !in.eof()) {
            getline(in,line);
            lines.push_back(line);
        }

    if(in.eof())
    {
        in.clear();
        in.seekg(0,ios_base::end);
    }

    auto current_line_no=int(lines.size())-1;

    string current_line;
    if(current_line_no>=0)
        current_line = lines[current_line_no];

    int position_within_current_line=0;
    if(in.tellg()>=0)
        position_within_current_line=current_line.length()-(int(in.tellg())-context_position);

    if (position_within_current_line<0)
        position_within_current_line=0;

    ostringstream lines_before;
    if(current_line_no==1) {
        lines_before << lines[0];
    } else if(current_line_no>1) {
        lines_before << lines[current_line_no-2] << endl << lines[current_line_no-1];
    }

    ostringstream lines_after;

    if(getline(in,line)) {
        lines_after << line;
    }
    if(getline(in,line)) {
        lines_after << endl << line;
    }

    in.seekg(context_position, ios_base::beg);
    string before=lines_before.str();
    string cur=current_line;
    string after=lines_after.str();

    ostringstream context;
    context << lines_before.str() << endl << current_line << endl;
    for(int i=0; i<position_within_current_line; ++i)
        context << " ";
    context << "^ - around here" << endl;
    context << lines_after.str();

    return ContextAroundPosition{lines_before.str(), current_line, lines_after.str(), position_within_current_line, current_line_no, context.str()};
}

string throw_parser_error(const string& filename, istream& in,const string& description) {
    ostringstream err_text;
    auto ctx = get_context(in);

    err_text << "Error parsing \"" << filename << "\" file" <<  " at line " << ctx.line_number << endl;

    err_text << description << ":";
    err_text << endl << endl << ctx.context << endl << endl;

    throw ParserError(err_text.str());
}

string throw_undefined_keyword(const string& filename, string keyword) {
    ostringstream err_text;
    err_text << "Error parsing \"" << filename << "\" file:" << endl;
    err_text << "keyword " << keyword << " is missing." << endl;

    throw ParserError(err_text.str());
}

string throw_error(const string& filename, string error) {
    ostringstream err_text;
    err_text << "Error parsing \"" << filename << "\" file:" << endl;
    err_text << error << endl;

    throw ParserError(err_text.str());
}

vector<string> split(const string& input_string, char where) {
    vector<string> res;

    std::string::size_type prev_pos = 0, pos = 0;
    while( (pos = input_string.find(where, pos)) != std::string::npos )
    {
        std::string substring( input_string.substr(prev_pos, pos-prev_pos) );

        if(substring.length()>0)
            res.push_back(substring);

        prev_pos = ++pos;
    }
    std::string substring( input_string.substr(prev_pos, pos-prev_pos) ); // Last word
    res.push_back(substring);
    return res;
}

string to_lower(const string& inp) {
    string res(inp);

    std::transform(res.begin(), res.end(), res.begin(), ::tolower);
    return res;
}

bool parse_bool(ifstream& in, string filename) {
    string user_input;
    in >> user_input;
    string user_input_l = to_lower(user_input);
    if(user_input_l == "true" || user_input_l == "1" || user_input_l == "yes")
        return true;
    else if(user_input_l == "false" || user_input_l == "0" || user_input_l == "no")
        return false;
    else
        throw_parser_error(filename, in, "Unexpected value\""+ user_input + "\"");

    return false;
}

InputParameters parse_input(const string& filename) {
    ifstream in(filename);

    if(!in)
        throw FileNotFound(filename);

    InputParameters par;

    string keyword;

    while (!in.eof()) {
        in >> keyword;
        if (in.eof())
            break;

        if (keyword[0] == '!' or keyword[0] == '#') {
            getline(in, keyword);
            continue;
        }

        //TODO: make a gentle warning when input parameter is .h5 and not text

//TODO: Complain if unit cell in hexagonal system, or in any other system really is not compatible
//with the angles in the input file
//not just average

        if (keyword == "INPUTS") {
        ////TODO: here remove bugs. Split names not by ' ' but by '\s'. trim string before doing so
        //possibly also by "some long pathname" and 'some long pathname', forget escapes for now
            string inputs;
            getline(in, inputs);
            par.input_files = split(inputs, ' ');
//            in >> par.data_filename_template;
        }else if (keyword == "OUTPUT")
            in >> par.output_name;
        else if (keyword == "SYMMETRY") {
            in >> par.symmetry;
            if(! isIn(par.symmetry, {"-1", "mmm", "m-3m", "1", "-3", "2/m:b", "2/m:c",
                                     "4/m", "4/mmm", "6/m", "6/mmm", "m-3", "-3:H", "-3:R", "-3m:R", "-3m:H"}))
                throw_parser_error(filename, in, "Unknown symmetry \""+ par.symmetry + "\"");
        }else if (keyword == "REJECT_OUTLIERS") {
            par.reject_outliers = parse_bool(in, filename);
        }else if(keyword == "SLICE"){
            par.slice = vector<float>(6,0);
            in >> par.slice[0] >> par.slice[1] >> par.slice[2] >> par.slice[3] >> par.slice[4] >> par.slice[5];
        }else if(keyword == "THRESHOLD")
            in >> par.threshold;
        else if(keyword == "ADD_CONSTANT")
            in >> par.add_constant;
        else if(keyword == "REPORT_PIXEL_RINT")
            par.report_pixel_rint = parse_bool(in, filename);
        else if(keyword == "REPORT_PIXEL_VARIANCE")
            par.report_pixel_variance = parse_bool(in, filename);
        else if(keyword == "REPORT_PIXEL_MULTIPLICITY")
            par.report_pixel_multiplicity = parse_bool(in, filename);
        else if(keyword == "SCALES") {
            //TODO: here remove bugs. Split names not by ' ' but by '\s'. trim string before doing so
            string inputs;
            getline(in, inputs);
            auto t = split(inputs, ' ');
            par.scales = vector<float>(t.size(), 0);
            transform(t.begin(), t.end(), par.scales.begin(), [](const string& in) {return ::atof(in.c_str());});
        }else if (keyword=="PUNCH_AND_FILL") {
            par.punch_and_fill = true;
            in >> par.r_punch >> par.r_fill;
        }else if(keyword=="BIN"){
            par.bin = true;
            in >> par.binning[0] >> par.binning[1] >> par.binning[2];
        }else if(keyword=="FFT"){
            par.fft = true;
        }else if(keyword=="CENTERING")
            in >> par.centering;
        else if(keyword=="TRIM_FOR_YELL")
            par.trim_for_yell = true;
        else {
            throw_parser_error(filename, in, "Unknown keyword \"" + keyword + "\"");
        }

        //Minor problem if file ends exactly with the last value, then this error fires
        if(in.fail() && in.eof()) {
            throw_parser_error(filename, in, "Unexpected end of file");
        }
        else if(in.fail()) {
            throw_parser_error(filename, in, "Incorrect arguments to keyword \"" + keyword + "\"");
        }
    }

    return par;
}
