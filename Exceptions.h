//
// Created by Arkadiy on 01/06/2017.
//

#ifndef MEERKAT2_AVERAGE_EXCEPTIONS_H
#define MEERKAT2_AVERAGE_EXCEPTIONS_H

#include <exception>
#include <string>

using std::string;

class FileNotFound : std::exception {
public:
    FileNotFound(string filename): filename(filename) {}
    ~FileNotFound() throw() {}
    string filename;
};

class ParserError : std::exception {
public:
    ParserError(string description): description(description) {}
    ~ParserError() throw() {}
    string description;
};

class UnknownFormat : std::exception {
public:
    UnknownFormat(){}
    ~UnknownFormat() throw() {}
};

class InvalidSlice : std::exception {
public:
    InvalidSlice(string description): description(description) {}
    ~InvalidSlice() throw() {}
    string description;
};

#endif //MEERKAT2_AVERAGE_EXCEPTIONS_H
