//
// Created by Arkadiy on 02/06/2017.
//

#ifndef MEERKAT2_AVERAGE_INTENSITYDATA_H
#define MEERKAT2_AVERAGE_INTENSITYDATA_H

#include <string>
#include <vector>
#include <assert.h>
#include <cmath>
#include <sstream>
#include <tools/kiss_fftnd.h>

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

struct Slice {
    bool everything;
    vector<float> limits;
    Slice(const vector<float> &limits):
            limits(limits),
            everything(false)
    {}

    Slice(): everything(true)
    {}
};

template<typename T>
inline T median(vector<T>& inp) {
    if(inp.size() == 0)
        return NAN;

    sort(inp.begin(), inp.end());
    return inp[inp.size()/2];
}

template<typename T>
inline T mean(vector<T>& inp) {
    T sum = 0;
    for(auto& v : inp)
        sum += v;

    return sum/inp.size();
}

inline int min(int a, size_t b) {
    return min(a, static_cast<int> (b));
}

//const Slice Slice::ALL

template<typename T>
class IntensityData {
public:
    IntensityData() {}

    static inline IntensityData<T> read(const string& filename, const Slice& slice) {
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

        //Read dataset and crystal data
        res.isDirect = readConstant<bool>(dataFile, "is_direct");
        res.unit_cell = readVector<double, 6>(dataFile, "unit_cell");

        res.lower_limits = readVector<double,3>(dataFile,"lower_limits");
        res.step_sizes = readVector<double, 3>(dataFile, "step_sizes");

        DataSet data = dataFile.openDataSet("data");
        hsize_t datasetDimesions[3];
        data.getSpace().getSimpleExtentDims( datasetDimesions, NULL);
        res.size = vector<size_t>(begin(datasetDimesions), end(datasetDimesions));

        //SLICE HERE
        if(slice.everything) {
            res.data = readVector<float>(dataFile, "data");
        }
        else{
            vector<int> start_pixels(3,0);
            vector<int> end_pixels(3,0);
            for(int i = 0; i<3; ++i){
                if(slice.limits[i] < res.lower_limits[i]-res.step_sizes[i]/2) {
                    stringstream err;

                    err << "lower limit error, requested: " << slice.limits[0] << ' ' << slice.limits[1] << ' ' << slice.limits[2];
                    err << " while dataset allows " << res.lower_limits[0] << ' ' << res.lower_limits[1] << ' ' << res.lower_limits[2];
                    throw InvalidSlice(err.str());
                }
                if(slice.limits[i + 3] > res.lower_limits[i]+res.step_sizes[i]*(res.size[i]-0.5)) {
                    stringstream err;

                    err << "upper limit error, requested: " << slice.limits[4] << ' ' << slice.limits[5] << ' ' << slice.limits[6];
                    err << " while dataset allows ";
                    for(int k = 0; k < 3; ++k)
                        err << res.lower_limits[k]+res.step_sizes[k]*res.size[k] << ' ';

                    throw InvalidSlice(err.str());
                }
                if(slice.limits[i] > slice.limits[i+3]) {
                    throw InvalidSlice("reverse slice is not allowed");
                }

                start_pixels[i] = round((slice.limits[i]-res.lower_limits[i])/res.step_sizes[i]);
                end_pixels[i]   = round((slice.limits[i+3]-res.lower_limits[i])/res.step_sizes[i]);
            }
            for(int i = 0; i < 3; ++i) {
                res.size[i] = end_pixels[i]-start_pixels[i]+1;
                res.lower_limits[i] += res.step_sizes[i]*start_pixels[i];
            }

            DataSet dataset = dataFile.openDataSet("data");
            DataSpace dataspace = dataset.getSpace();
            vector<hsize_t> memspaceOffset(start_pixels.begin(), start_pixels.end());
            vector<hsize_t> spaceCount(res.size.begin(), res.size.end());
            dataspace.selectHyperslab(H5S_SELECT_SET, spaceCount.data(), memspaceOffset.data());

            res.data = vector<T>(res.size[0]*res.size[1]*res.size[2]);

            DataSpace memspace(3, spaceCount.data());
            dataset.read(res.data.data(), getH5Type<T>(), memspace, dataspace);
        }

//    metricTensor = metricTensorFromUnitCell(unitCell, !isDirect);
//        DSetCreatPropList cparms = data.getCreatePlist();
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

    template<typename T2>
    static inline IntensityData empty(const IntensityData<T2>& inp) {
        IntensityData res;
        res.size = inp.size;
        res.isDirect = inp.isDirect;
        res.lower_limits = inp.lower_limits;
        res.step_sizes = inp.step_sizes;
        res.unit_cell = inp.unit_cell;
        res.data = vector<T>(inp.data.size(), 0);
        return res;
    }

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

    void scale(double sc) {
        for(auto& d : data)
            d *= sc;
    }

    void accumulate(IntensityData& inp, double scale = 1) {
        for(auto d1 = data.begin(), d2 = inp.data.begin(); d1 != data.end(); ++d1, ++d2)
        {
            *d1 += (*d2)*scale;
        }
    }

    void punch_and_fill(double r_punch, double r_fill, const string& centering) {
        //TODO: deal with the cases of non-cubic systems. For instance in hexagonal systems the radius will be dependent on the metric tensor.
        //TODO: also figure out what to do in cases when pixels are different in different directions (three radii??).
        //Will need an nd array class with slices and strides like in numpy
        //Will need slices for Intensity Data which take care of all the lower limits and step sizes
        //Will need distances in the metric space
        //TODO: this code is horrific. It is not a performance bottleneck, so can be improved if virtual slice and functions for dist are introduced.
        double eps = 1e-9;
        auto ul = upper_limits();

        int r = round(r_fill + 1);
        int sz = 2*r+1;
        vector<T> backgdound_pixels(sz*sz*sz);

        for(int h = ceil(lower_limits[0])-1; h < ul[0]+1; ++h)
            for(int k = ceil(lower_limits[1])-1; k < ul[1]+1; ++k)
                for(int l = ceil(lower_limits[2])-1; l < ul[2]+1; ++l) {
                    double hi_c = (h-lower_limits[0])/step_sizes[0];
                    double ki_c = (k-lower_limits[1])/step_sizes[1];
                    double li_c = (l-lower_limits[2])/step_sizes[2];

                    if(centering == "F") {
                        if(! ( abs(h%2)==abs(k%2) and abs(h%2)==abs(l%2) ))
                            continue;
                    }

                    backgdound_pixels.resize(0);
                    for(int dhi = -r, hi = round(hi_c)+dhi; dhi < r; ++dhi, ++hi)
                        if(hi >= 0 && hi < size[0])
                            for(int dki = -r, ki = round(ki_c)+dki; dki < r; ++dki, ++ki)
                                if(ki >= 0 && ki < size[1])
                                    for(int dli = -r, li = round(li_c)+dli; dli < r; ++dli, ++li)
                                        if(li >= 0 && li < size[2]) {
                                            double r = sqrt(dhi*dhi+dki*dki+dli*dli);
                                            if(r <= r_fill && r >= r_punch) {
                                                auto val = at(hi, ki, li);
                                                if(! isnan(val))
                                                    backgdound_pixels.push_back(val);
                                            }
                                        }

                    auto fill_val = median(backgdound_pixels);

                    for(int dhi = -r, hi = round(hi_c)+dhi; dhi < r; ++dhi, ++hi)
                        if(hi >= 0 && hi < size[0])
                            for(int dki = -r, ki = round(ki_c)+dki; dki < r; ++dki, ++ki)
                                if(ki >= 0 && ki < size[1])
                                    for(int dli = -r, li = round(li_c)+dli; dli < r; ++dli, ++li)
                                        if(li >= 0 && li < size[2]) {
                                            double r = sqrt(dhi*dhi+dki*dki+dli*dli);
                                            if(r < r_punch) {
                                                at(hi, ki, li) = fill_val;
                                            }
                                        }
                }
    }

    vector<double> upper_limits() {
        vector<double> res(3,0);
        for(int i = 0; i < 3; ++i) {
            res[i] = size[i]*step_sizes[i]+lower_limits[i];
        }
        return res;
    }

    vector<int> cenre_ind() {
        vector<int> res(3);
        for(int i = 0; i < 3; ++i)
            res[i] = round(-lower_limits[i]/step_sizes[i]);

        return res;
    }

    int ind2ind(int hi, int ki, int li) {
        return (hi*size[1]+ki)*size[2]+li;
    }

    T& at(int hi, int ki, int li) {
        return data[ind2ind(hi, ki, li)];
    }

    IntensityData<T> binned(vector<int> binning) {
        vector<int> starting_indices(3);
        vector<double> new_lower_limits(3);
        vector<double> new_step_sizes(3);
        vector<size_t> new_size(3);
        for(int i = 0; i < 3; ++i) {
            assert(binning[i] % 2 == 1);

            int central_pixel_index = round(-lower_limits[i]/step_sizes[i]);
            int available_pixels = central_pixel_index + binning[i]/2 - 1;
            int starting_index = central_pixel_index - round(ceil(static_cast<double> (available_pixels)/binning[i])*binning[i]);
            starting_indices[i] = starting_index;

            new_lower_limits[i] = lower_limits[i]+starting_indices[i]*step_sizes[i];

            int effective_pixels = size[i]-starting_indices[i];
            new_size[i] = ceil(static_cast<double> (effective_pixels)/binning[i]);
            new_step_sizes[i] = step_sizes[i]*binning[i];
        }

        IntensityData<T> res;
        res.size = new_size;
        res.lower_limits = new_lower_limits;
        res.step_sizes = new_step_sizes;
        res.isDirect = isDirect;
        res.unit_cell = unit_cell;
        res.data.resize(new_size[0]*new_size[1]*new_size[2]);

        for(int hr = 0, hic = starting_indices[0]; hr < new_size[0]; ++hr, hic+=binning[0])
            for(int kr = 0, kic = starting_indices[1]; kr < new_size[1]; ++kr, kic+=binning[1])
                for(int lr = 0, lic = starting_indices[2]; lr < new_size[2]; ++lr, lic+=binning[2]) {
                    T accum = 0;
                    int Npix = 0;

                    for(int hi = max(hic-binning[0]/2,0); hi < min(hic+binning[0]/2+1, size[0]); ++hi)
                        for(int ki = max(kic-binning[1]/2,0); ki < min(kic+binning[1]/2+1, size[1]); ++ki)
                            for(int li = max(lic-binning[2]/2,0); li < min(lic+binning[2]/2+1, size[2]); ++li) {
                                accum += at(hi,ki,li);
                                Npix += 1;
                            }
                    res.at(hr, kr, lr) += accum/Npix;
                }

        return res;
    }

    void FFT() {
        int dims[] = {0,0,0};
        auto trim = [](int dim) {
            if(dim == 1)
                return dim;
            else
                return dim - (dim % 2);
        };

        for(int i=0; i<3; ++i)
            dims[i] = trim(size[i]);

        kiss_fftnd_cfg  cfg = kiss_fftnd_alloc(dims, 3, 0, NULL, NULL);
        auto new_npix = dims[0]*dims[1]*dims[2];
        vector<kiss_fft_cpx> buffer(new_npix, {0,0});

        int bi = 0;
        vector<int> signs(3);
        auto cen = cenre_ind();
        transform(cen.begin(), cen.end(), signs.begin(), [](int n){
            if(n % 2 == 0)
                return 1;
            else
                return -1;
        });

        for(int hi = 0, hsign = signs[0]; hi < dims[0]; ++hi, hsign*=-1)
            for(int ki = 0, ksign = signs[1]; ki < dims[1]; ++ki, ksign*=-1)
                for(int li = 0, lsign = signs[2]; li < dims[2]; ++li, lsign*=-1) {
                    auto val = at(hi, ki, li);
                    if(!isnan(val))
                        buffer[bi].r = val*hsign*ksign*lsign;

                    ++bi;
                }

        kiss_fftnd(cfg, buffer.data(), buffer.data());

        data.resize(new_npix);

        for(int i = 0; i < 3; ++i) {
            size[i] = dims[i];
            auto t = lower_limits[i];
            lower_limits[i] = -0.5/step_sizes[i];
            step_sizes[i] = -0.5/t;
            isDirect = !isDirect;
        }

        for(int hi = 0, hsign = signs[0]; hi < dims[0]; ++hi, hsign*=-1)
            for(int ki = 0, ksign = signs[1]; ki < dims[1]; ++ki, ksign*=-1)
                for(int li = 0, lsign = signs[2]; li < dims[2]; ++li, lsign*=-1) {
                    at(hi,ki,li) = buffer[ind2ind(hi,ki,li)].r*hsign*ksign*lsign;
                }
    }
private:
};



#include "Exceptions.h"





#endif //MEERKAT2_AVERAGE_INTENSITYDATA_H
