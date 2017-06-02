//
// Created by Arkadiy on 02/06/2017.
//

#include "IntensityData.h"

bool isYellFormat(H5File& file) {
    std::string format;
    try {
        StrType datatype(0, H5T_VARIABLE);
        DataSpace dataspace(H5S_SCALAR);
        DataSet datset = file.openDataSet("format");

        datset.read(format, datatype, dataspace);
    } catch(H5::FileIException) {
        return false;
    }

    if(format!="Yell 1.0")
        return false;

    else
        return true;
}

void writeYellFormatString(H5File& file) {
    string format = "Yell 1.0";

    H5::StrType h5stringType(H5::PredType::C_S1, H5T_VARIABLE); // + 1 for trailing zero
    H5::DataSet ds = file.createDataSet("format", h5stringType, H5::DataSpace(H5S_SCALAR));
    ds.write(format, h5stringType);
}