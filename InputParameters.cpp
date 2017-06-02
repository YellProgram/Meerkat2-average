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

bool reached_eof(const istream& in) {
    return in.rdstate() & istream::eofbit;
}

bool failed(const istream& in) {
    return in.rdstate() & (istream::badbit | istream::failbit);
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


InputParameters parse_input(const string& filename) {
    ifstream in(filename);

    if(!in)
        throw FileNotFound(filename);

    InputParameters par;

    string keyword;

    while (!in.eof()) {
        in >> keyword;
        if (reached_eof(in))
            break;

        if (keyword[0] == '!' or keyword[0] == '#') {
            getline(in, keyword);
            continue;
        }

        if (keyword == "INPUTS") {
            string inputs;
            getline(in, inputs);
            par.input_files = split(inputs, ' ');
//            in >> par.data_filename_template;
        }else if (keyword == "OUTPUT")
            in >> par.output_name;
        else if (keyword == "SYMMTERY") {
            in >> par.symmetry;
            if(! isIn(par.symmetry, {"-1", "mmm", "m3m"}))
                throw_parser_error(filename, in, "Unknown symmetry\""+ par.symmetry + "\"");
        }

        else if (keyword == "REJECT_OUTLIERS") {
            string user_input;
            in >> user_input;
            string user_input_l = to_lower(user_input);
            if(user_input_l == "true" || user_input_l == "1" || user_input_l == "yes")
                par.reject_outliers = true;
            else if(user_input_l == "false" || user_input_l == "0" || user_input_l == "no")
                par.reject_outliers = false;
            else
                throw_parser_error(filename, in, "Unexpected value\""+ user_input + "\"");
        }
        else if(keyword == "THRESHOLD")
            in >> par.threshold;
        else {
            throw_parser_error(filename, in, "Unknown keyword \"" + keyword + "\"");
        }

        //minor problem if file ends exactly with the last value, then this error fires
        if(reached_eof(in)) {
            throw_parser_error(filename, in, "Unexpected end of file");
        }
        else if(failed(in)) {
            throw_parser_error(filename, in, "Incorrect arguments to keyword \"" + keyword + "\"");
        }

    }

    return par;
}