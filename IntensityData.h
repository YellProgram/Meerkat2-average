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
#include <map>

#include "Exceptions.h"
#include <H5Cpp.h>

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
vector<vector<A> > readMatrix(H5File& f, const string& datasetName) {
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
vector<A> readVector(H5File& f, const string& datasetName) {
    DataSet dataset = f.openDataSet(datasetName);

    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    assert(rank==1);
    hsize_t count[1] = {datasetSize};
    dataspace.getSimpleExtentDims(count);
    assert(count[0]==datasetSize);

    A out[datasetSize];
    dataset.read( out, getH5Type<A>());
    vector<A>result (begin(out), end(out));
    return result;
}

template<typename A>
A readConstant(H5File& f, const string& datasetName) {
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
    IntensityData():
            Rint(NAN)
    {}

    static inline IntensityData<T> read(const string& filename, const Slice& slice) {
        IntensityData<T> res;

        auto readProps = FileAccPropList::DEFAULT;
        int mdc_nelmts;
        size_t rdcc_nelmts, rdcc_nbytes;
        double rdcc_w0;
        readProps.getCache(mdc_nelmts, rdcc_nelmts, rdcc_nbytes, rdcc_w0);
        readProps.setCache(mdc_nelmts, 269251, 1024*1024*500, rdcc_w0); //Magic numbers to tune performance
        H5File dataFile = H5File(filename, H5F_ACC_RDONLY, FileCreatPropList::DEFAULT, readProps);

        if(!isYellFormat(dataFile))
            throw UnknownFormat();

        //Read dataset and crystal data
        res.isDirect = readConstant<bool>(dataFile, "is_direct");
        res.unit_cell = readVector<double, 6>(dataFile, "unit_cell");

        res.lower_limits = readVector<double, 3>(dataFile,"lower_limits"); //TODO: fix bug with returns here. They do not run with -O2 or -O3 flags

        res.step_sizes = readVector<double, 3>(dataFile, "step_sizes");

        DataSet data = dataFile.openDataSet("data");
        hsize_t datasetDimesions[3];
        data.getSpace().getSimpleExtentDims( datasetDimesions, NULL);
        res.size = vector<size_t>(begin(datasetDimesions), end(datasetDimesions));

        //SLICE HERE
        if(slice.everything) {
            res.data = readVector<T>(dataFile, "data");
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

        
        assert(res.lower_limits.size()==3); //for debug purposes. TODO: remove after teh bug with lower_limits is fixed
        return res;
    }

    template<typename T2>
    void saveData(string filename, vector<T2>& dataset) {
        H5File file( filename, H5F_ACC_TRUNC );

        creadeAndWriteDataset<float>(file, "data", dataset, size);
        creadeAndWriteDataset(file, "lower_limits", lower_limits);
        writeConstant(file, "is_direct", isDirect);
        if(!isnan(Rint))
            writeConstant(file, "Rint", Rint);

        creadeAndWriteDataset(file, "step_sizes", step_sizes);
        creadeAndWriteDataset(file, "unit_cell", unit_cell);
        writeYellFormatString(file);
    }

    void save(string filename) {
        saveData(filename, data);

        for(auto& it : other_datasets) {
            auto pos_extension = filename.find(".");
            auto prefix = filename.substr(0,pos_extension);
            auto out_filename = prefix + "_" + it.first + ".h5";
            saveData(out_filename, it.second);
        }
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
    map<string, vector<float>> other_datasets;
    float Rint;

    void scale(double sc) {
        for(auto& d : data)
            d *= sc;
    }

    void set_nans_to_zeros() {
        for(auto d1 = data.begin(); d1 != data.end(); ++d1)
            if (isnan(*d1))
                *d1=0;
    }

    void divide_by(IntensityData& inp) {
        for(auto d1 = data.begin(), d2 = inp.data.begin(); d1 != data.end(); ++d1, ++d2)
        {
            *d1 /= *d2;
        }
    }

    void accumulate(IntensityData& inp, double scale = 1, bool accumulate_nans = true) {
        for(auto d1 = data.begin(), d2 = inp.data.begin(); d1 != data.end(); ++d1, ++d2)
        {
            auto val = *d2;
            if(accumulate_nans || !isnan(val)) //Accumulate propagating nans if the flag tells so, otherwise only add if the value is non nan.
                *d1 += val*scale;
        }
    }

    ///Add the value of scale to each pixel which is measured (not equal to nan) in the dataset.
    void accumulate_non_nan_scales(IntensityData& inp, double scale = 1) {
        for(auto d1 = data.begin(), d2 = inp.data.begin(); d1 != data.end(); ++d1, ++d2)
        {
            if(!isnan(*d2))
                *d1 += scale;
        }
    }

    void punch_and_fill(double r_punch, double r_fill, const string& centering) {
        //TODO: deal with the cases of non-cubic systems. For instance in hexagonal systems the radius will be dependent on the metric tensor.
        //TODO: also figure out what to do in cases when pixels are different in different directions (three radii??).
        //Will need an nd array class with slices and strides like in numpy
        //Will need slices for Intensity Data which take care of all the lower limits and step sizes
        //Will need distances in the metric space
        //TODO: this code is not very beautiful. It is not a performance bottleneck, so can be improved if virtual slice and functions for dist are introduced.

        //TODO: fix a bug do not fill in regions which were nans
        double eps = 1e-9;
        auto ul = upper_limits();

        int r = round(max(r_fill,r_punch) + 1);
        int sz = 2*r+1;
        vector<T> backgdound_pixels(sz*sz*sz);

        for(int h = ceil(lower_limits[0])-1; h < ul[0]+1; ++h)
            for(int k = ceil(lower_limits[1])-1; k < ul[1]+1; ++k)
                for(int l = ceil(lower_limits[2])-1; l < ul[2]+1; ++l) {
                    double hi_c = (h-lower_limits[0])/step_sizes[0];
                    double ki_c = (k-lower_limits[1])/step_sizes[1];
                    double li_c = (l-lower_limits[2])/step_sizes[2];

                    if(centering == "F") {
                        if(! ( abs(h%2)==abs(k%2) && abs(h%2)==abs(l%2) ))
                            continue;
                    }else if(centering == "I") {
                        if(! ((h+k+l)%2 == 0))
                            continue;
                    }else if(centering == "R") {
                        if(! ((-h+k+l)%3 == 0))
                            continue;
                    }else if(centering == "C") {
                      if(! ((h+k) % 2 == 0))
                        continue;
                    }else if(centering == "B") {
                        if(! ((h+l) % 2 == 0))
                            continue;
                    }else if(centering == "A") {
                        if(! ((k+l) % 2 == 0))
                            continue;
                    } //TODO: add checks for unknown centerings. Need R reverse and obverse? probably not

                    backgdound_pixels.resize(0);
                    for(int dhi = -r, hi = round(hi_c)+dhi; dhi <= r; ++dhi, ++hi)
                        if(hi >= 0 && hi < size[0])
                            for(int dki = -r, ki = round(ki_c)+dki; dki <= r; ++dki, ++ki)
                                if(ki >= 0 && ki < size[1])
                                    for(int dli = -r, li = round(li_c)+dli; dli <= r; ++dli, ++li)
                                        if(li >= 0 && li < size[2]) {
                                            double r = sqrt(dhi*dhi+dki*dki+dli*dli);
                                            if(r <= r_fill && r >= r_punch) {
                                                auto val = at(hi, ki, li);
                                                if(! isnan(val))
                                                    backgdound_pixels.push_back(val);
                                            }
                                        }

                    auto fill_val = median(backgdound_pixels);

                    for(int dhi = -r, hi = round(hi_c)+dhi; dhi <= r; ++dhi, ++hi)
                        if(hi >= 0 && hi < size[0])
                            for(int dki = -r, ki = round(ki_c)+dki; dki <= r; ++dki, ++ki)
                                if(ki >= 0 && ki < size[1])
                                    for(int dli = -r, li = round(li_c)+dli; dli <= r; ++dli, ++li)
                                        if(li >= 0 && li < size[2]) {
                                            double r = sqrt(dhi*dhi+dki*dki+dli*dli); //TODO: replace this with proper distance in reciprocal angstroems
                                            if(r < r_punch && !isnan(at(hi, ki, li))) {
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
        // Must use odd binning to maintain a clear center point.
        assert(binning[i] % 2 == 1);

        // Calculate the number of new bins (includes partial bins).
        size_t num_bins = static_cast<size_t>(ceil(static_cast<double>(size[i]) / binning[i]));
        new_size[i] = num_bins;

        // Calculate the total number of pixels needed to form these 'num_bins' exactly.
        size_t needed_pixels = num_bins * binning[i];

        // Determine the Starting Index (SI).
        // SI is the number of 'virtual' pixels added to the negative side
        // to ensure center-alignment (e.g., -2 for 601 pixels and binning 5).
        int added_pixels = needed_pixels - size[i];
        int starting_index = -added_pixels / 2;
        starting_indices[i] = starting_index;

        // The new physical lower limit is shifted by the 'virtual' starting pixels.
        new_lower_limits[i] = lower_limits[i] + starting_indices[i] * step_sizes[i];

        new_step_sizes[i] = step_sizes[i] * binning[i];
    }

    IntensityData<T> res;
    res.size = new_size;
    res.lower_limits = new_lower_limits;
    res.step_sizes = new_step_sizes;
    res.isDirect = isDirect;
    res.unit_cell = unit_cell;
    res.data.resize(new_size[0] * new_size[1] * new_size[2]);

    // Outer loops iterate over the new, smaller size (hr, kr, lr).
    for(int hr = 0; hr < new_size[0]; ++hr)
        for(int kr = 0; kr < new_size[1]; ++kr)
            for(int lr = 0; lr < new_size[2]; ++lr) {

                // Calculate the block's START and END indices in the original *virtual* space.
                int h_virtual_start = starting_indices[0] + hr * binning[0];
                int k_virtual_start = starting_indices[1] + kr * binning[1];
                int l_virtual_start = starting_indices[2] + lr * binning[2];

                int h_virtual_end = h_virtual_start + binning[0];
                int k_virtual_end = k_virtual_start + binning[1];
                int l_virtual_end = l_virtual_start + binning[2];

                T accum = 0;
                int Npix = 0;

                // Inner loops iterate only over the part of the bin that exists
                // within the original array's bounds [0, size[i]).
                for(int hi = max(h_virtual_start, 0); hi < min(h_virtual_end, (int)size[0]); ++hi)
                    for(int ki = max(k_virtual_start, 0); ki < min(k_virtual_end, (int)size[1]); ++ki)
                        for(int li = max(l_virtual_start, 0); li < min(l_virtual_end, (int)size[2]); ++li) {
                            accum += at(hi,ki,li);
                            Npix += 1;
                        }

                if (Npix > 0) {
                    res.at(hr, kr, lr) = accum / Npix;
                } else {
                    res.at(hr, kr, lr) = 0;
                }
            }

    return res;
}

    void prepare_for_fft() {
        //cut last rows

        vector<size_t> dims(3);

        auto trim = [](int dim) {
            if(dim == 1)
                return dim;
            else
                return dim - (dim % 2);
        };

        for(int i=0; i<3; ++i)
            dims[i] = trim(size[i]);

        //trim last rows
        int bi = 0;
        for(int hi = 0; hi < dims[0]; ++hi)
            for(int ki = 0; ki < dims[1]; ++ki)
                for(int li = 0; li < dims[2]; ++li) {
                    data[bi] = at(hi, ki, li);

                    ++bi;
                }

        size = dims;
        data.resize(size[0]*size[1]*size[2]);
    }

    void set_nans_to_zero() {
        for(int bi = 0; bi < data.size(); ++bi) {
                    auto val = data[bi];
                    if(!isnan(val))
                        data[bi] = val;
                    else
                        data[bi] = 0;

                    ++bi;
                }

    }

    void add_constant(double c) {
        for(auto& v : data)
            v+=c;
    }

    void FFT() {
        vector<int> dims(size.begin(), size.end());

        kiss_fftnd_cfg  cfg = kiss_fftnd_alloc(dims.data(), 3, 0, NULL, NULL);
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
                    buffer[bi].r = at(hi, ki, li)*hsign*ksign*lsign;
                    ++bi;
                }

        kiss_fftnd(cfg, buffer.data(), buffer.data());

        for(int i = 0; i < 3; ++i) {
            if (dims[i]==1) {
                lower_limits[i]=0;
                step_sizes[i]=1;
            } else {
                auto const t = lower_limits[i];
                lower_limits[i] = -0.5/step_sizes[i];
                step_sizes[i] = -0.5/t;
            }

        }
        isDirect = !isDirect;

        for(int hi = 0, hsign = signs[0]; hi < dims[0]; ++hi, hsign*=-1)
            for(int ki = 0, ksign = signs[1]; ki < dims[1]; ++ki, ksign*=-1)
                for(int li = 0, lsign = signs[2]; li < dims[2]; ++li, lsign*=-1) {
                    at(hi,ki,li) = buffer[ind2ind(hi,ki,li)].r*hsign*ksign*lsign;
                }
    }
private:
};



#include "Exceptions.h"

//TODO: add all possible centering types (C, I)



#endif //MEERKAT2_AVERAGE_INTENSITYDATA_H
