//
// Created by Arkadiy on 02/06/2017.
//

#ifndef MEERKAT2_AVERAGE_INTENSITYDATA_H
#define MEERKAT2_AVERAGE_INTENSITYDATA_H

#include <string>
#include <vector>
#include <assert.h>

#include "Exceptions.h"
#include "H5Cpp.h"

using namespace std;
using namespace H5;



template<typename T>
inline DataType getH5Type() {}

template<>
inline DataType getH5Type<int> () {
    return PredType::NATIVE_INT;
}
template<>
inline DataType getH5Type<bool> () {
    return PredType::NATIVE_HBOOL;
}
template<>
inline DataType getH5Type<double> () {
    return PredType::NATIVE_DOUBLE;
}
template<>
inline DataType getH5Type<float> () {
    return PredType::NATIVE_FLOAT;
}

template<typename A, hsize_t Nx, hsize_t Ny>
vector<vector<A> > readMatrix(H5File f, const string& datasetName) {
    DataSet dataset = f.openDataSet(datasetName);

    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    assert(rank==2);
    hsize_t count[2] = {Nx, Ny};
    dataspace.getSimpleExtentDims(count);
    assert(count[0]==Nx && count[1]==Ny);

    A out[Nx][Ny];
    dataset.read( out, getH5Type<A>());

    vector<vector<A> > res(Nx);
    for(int i=0; i<Nx; ++i)
        res[i]=vector<A> (begin(out[i]),end(out[i]));

    return res;
}

template<typename A, hsize_t datasetSize>
vector<A> readVector(H5File f, const string& datasetName) {
    DataSet dataset = f.openDataSet(datasetName);

    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    assert(rank==1);
    hsize_t count[1] = {datasetSize};
    dataspace.getSimpleExtentDims(count);
    assert(count[0]==datasetSize);

    A out[datasetSize];
    dataset.read( out, getH5Type<A>());

    return vector<A>(begin(out),end(out));
}

template<typename A>
A readConstant(H5File f, const string& datasetName) {
    //Due to (probably) bug, dataset.read overwrites two values for a scalar
    //don't have time to trace it
    A res[2];

    DataSet dataset = f.openDataSet(datasetName);
    DataSpace dataspace = dataset.getSpace();
    assert(dataspace.getSimpleExtentNdims()==0);

    dataset.read(res, getH5Type<A>(),DataSpace(H5S_SCALAR),DataSpace(H5S_SCALAR));

    return res[0];
}

bool isYellFormat(H5File& file);

template<typename A>
vector<A> readVector(H5File& f, const string& datasetName) {
    DataSet dataset = f.openDataSet(datasetName);

    DataSpace dataspace = dataset.getSpace();

    hsize_t datasetSize = dataspace.getSimpleExtentNpoints();

    std::vector<A> out(datasetSize);
    dataset.read( out.data(), getH5Type<A>());

    return out;
}

template <typename T>
void creadeAndWriteDataset(H5File& file, string datasetName, vector<T> data, vector<size_t> inp_dims) {

    vector<hsize_t> dims(inp_dims.begin(), inp_dims.end());
    DataSpace dataspace( dims.size(), dims.data() );
    DataSet dataset = file.createDataSet( datasetName, getH5Type<T>(), dataspace );
    dataset.write( data.data(), getH5Type<T>() );
}

template <typename T>
void creadeAndWriteDataset(H5File& file, string datasetName, vector<T> data) {
    creadeAndWriteDataset(file, datasetName, data, {data.size()});
}

template <typename T>
void writeConstant(H5File& file, string datasetName, T data) {
    H5::DataSet ds = file.createDataSet(datasetName, getH5Type<T>(), H5::DataSpace(H5S_SCALAR));
    ds.write(&data, getH5Type<T>());
}


void writeYellFormatString(H5File& file);

//vector<vector<double>> metricTensorFromUnitCell(
//        const vector<double>& unitCell,
//        bool invert)
//{
//    double a = unitCell[0];
//    double b = unitCell[1];
//    double c = unitCell[2];
//    double alpha = deg2rad(unitCell[3]);
//    double beta  = deg2rad(unitCell[4]);
//    double gamma = deg2rad(unitCell[5]);
//    QMatrix4x4 t (a*a, a*b*cos(gamma), a*c*cos(beta),0,
//                  b*a*cos(gamma),b*b, b*c*cos(alpha),0,
//                  c*a*cos(beta),c*b*cos(alpha),c*c,0,
//                  0,0,0,1);
//    if(invert)
//        t=t.inverted();
//
//    auto tt=t.data();
//    return {{tt[0],tt[1],tt[2]},
//            {tt[4],tt[5],tt[6]},
//            {tt[8],tt[9],tt[10]}};
//}

// A class to interface with Yell files in hdf5 format
// should be flexible enough to be used here, in DensityViewer and in Meerkat
// the point of flexibility is in having two options - in memory data processing
// and on the hard drive for memory economy

// for this script it should allow to read the data from hard drive
// create one from scratch, save it to hard drive

// conveniently create new empty one from example
// LATER: subsection
// LATER: fft (if I need FFT i will need temporarily create complex double dataset

template<typename T>
class IntensityData {
public:
    IntensityData() {}

    static inline IntensityData<T> read(string filename) {
        IntensityData<T> res;

        auto readProps = FileAccPropList::DEFAULT;
        int mdc_nelmts;
        size_t rdcc_nelmts, rdcc_nbytes;
        double rdcc_w0;
        readProps.getCache(mdc_nelmts, rdcc_nelmts, rdcc_nbytes, rdcc_w0);
        readProps.setCache(mdc_nelmts, 269251, 1024*1024*500, rdcc_w0); //Magic numbers to tune performance
        H5File dataFile = H5File(filename, H5F_ACC_RDONLY, FileCreatPropList::DEFAULT, readProps);

        if(not isYellFormat(dataFile))
            throw UnknownFormat();

        DataSet data = dataFile.openDataSet("data");

        //Read dataset and crystal data
        res.isDirect = readConstant<bool>(dataFile, "is_direct");

        res.lower_limits = readVector<double,3>(dataFile,"lower_limits");
        res.step_sizes = readVector<double, 3>(dataFile, "step_sizes");
        res.unit_cell = readVector<double, 6>(dataFile, "unit_cell");
//    metricTensor = metricTensorFromUnitCell(unitCell, !isDirect);

        DSetCreatPropList cparms = data.getCreatePlist();

        hsize_t datasetDimesions[3];
        data.getSpace().getSimpleExtentDims( datasetDimesions, NULL);

        res.size = vector<size_t>(begin(datasetDimesions), end(datasetDimesions));

        res.data = readVector<float>(dataFile, "data");
        return res;
    }

     void save(string filename) {
        H5File file( filename, H5F_ACC_TRUNC );

        creadeAndWriteDataset<float>(file, "data", data, size);
        creadeAndWriteDataset(file, "lower_limits", lower_limits);
        writeConstant(file, "is_direct", false);
        creadeAndWriteDataset(file, "step_sizes", step_sizes);
        creadeAndWriteDataset(file, "unit_cell", unit_cell);
        writeYellFormatString(file);
    }


    static IntensityData empty(const IntensityData& inp);


    vector<size_t> size;
//    bool in_memory;

    vector<T> data;
//double at(int x, int y, int z) {return data[(x*size[1]+y)*size[2]+z];}
//vector<double> ind2hkl(const vector<int> & indices);
//    double lowerLimit(int i);
//    double upperLimit(int i);
//    double stepSize(int i);
    bool isDirect;
//    H5File dataFile;
//    DataSet data;

    vector<double> lower_limits;
    vector<double> step_sizes;
//    vector<vector<double> > metric_tensor;
    vector<double> unit_cell;
private:
};



#include "Exceptions.h"





#endif //MEERKAT2_AVERAGE_INTENSITYDATA_H
