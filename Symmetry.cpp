//
// Created by Arkadiy on 01/06/2017.
//

#include <iostream>
#include "Symmetry.h"
#include "InputParameters.h"


vector<Matrix3i> expand_generators(const vector<Matrix3i>& generators) {

    auto I = Matrix3i::Identity();
    vector<Matrix3i> res = {I};

    for(auto& g : generators) {
        auto last_iteration = res;
        auto gn = g;

        while(gn != I) {
            for(auto& h : last_iteration)
                res.push_back(gn*h);

            gn = gn*g;
        }
    }

    auto last = unique(res.begin(), res.end());
    res.erase(last, res.end());
    return res;
}

Matrix3i el(Vector3i x, Vector3i y, Vector3i z) {
    Matrix3i res;
    res << x, y, z;
    return res;
}

vector<Matrix3i> expand_symmetry(const string& symmetry) {
    vector<Matrix3i> res;

    Vector3i x, y, z;
    x << 1,0,0;
    y << 0,1,0;
    z << 0,0,1;

    if (symmetry=="m3m") {
        return expand_generators({el(-x, y, z),
                                  el( x,-y, z),
                                  el( x, y,-z),
                                  el( y, z, x),
                                  el( y, x, z),
                                  el(-y, x, z)});
    }else if(symmetry == "mmm") {
        return expand_generators({el(-x, y, z),
                                  el( x,-y, z),
                                  el( x, y,-z)});
    }else if(symmetry == "-1") {
        return expand_generators({el(-x, -y, -z)});
    }else if(symmetry == "1") {
        return expand_generators({el(x, y, z)});
    }else {
        assert(false);
    }
    return res;
}

inline size_t ind2ind(const Vector3i& r, const vector<size_t>& size,const vector<int>& centre) {
    return ((r[0]+centre[0])*size[1]+(r[1]+centre[1]))*size[2]+(r[2]+centre[2]);
}

inline float mean(const vector<float>::iterator &start, const vector<float>::iterator &end) {
    float average = 0;
    for_each(start,
             end,
             [&](const float& i){
                 average += i;
             });
    average = average/distance(start, end);
    return average;
}

inline float median(const vector<float>::iterator &start, const vector<float>::iterator &end) {
    sort(start, end);
    auto sz = distance(start, end);
    return *(start+sz/2);
}

inline const vector<float>::iterator reject_outliers(const vector<float>::iterator & start, const vector<float>::iterator & end, float threshold) {
    auto sz = distance(start, end);
    float med = median(start, end);

    vector<float> differences(sz);
    transform(start, end, differences.begin(), [&](const float& i) {return i - med;} );

    float median_distance = median(differences.begin(), differences.end());

    return remove_if(start, end, [&](const float& i){return abs(i - med)>median_distance*threshold;});
}

inline float average_set_of_reflections(const vector<float>::iterator &equivalent_intensities,  vector<float>::iterator &intensities_last, const InputParameters& par) {
    if(par.reject_outliers && distance(equivalent_intensities, intensities_last)>2) {
        intensities_last = reject_outliers(equivalent_intensities, intensities_last, par.threshold);
    }

    return mean(equivalent_intensities, intensities_last);
}

void average(IntensityData<float>& inp, IntensityData<float>& res, const InputParameters& par) {
    vector<int> centre(3);
    for(int i=0; i<3; ++i) {
        centre[i] = round(-inp.lower_limits[i]/inp.step_sizes[i]);
    }

    auto is_reconstructed = IntensityData<bool>::empty(inp);

    auto size = inp.size;

    //TODO: check grid is ok with symmetry, or otherwise ignore the portion which is outside when reconstructing
    //TODO: rewrite this loop in a linear fashion??
    vector<Matrix3i> symmetry_elements = expand_symmetry(par.symmetry);
    vector<size_t> equivalent_indices(symmetry_elements.size());
    vector<float> equivalent_intensities(symmetry_elements.size());
    Vector3i r;
    for(r(0) = -centre[0]; r(0) < (long)(size[0])-centre[0]; ++r(0))
        for(r(1) = -centre[1]; r(1) < (long)(size[1])-centre[1]; ++r(1))
            for(r(2) = -centre[2]; r(2) < (long)(size[2])-centre[2]; ++r(2))
                if(!is_reconstructed.data[ind2ind(r, size, centre)]) {
                    //select unique indices
                    transform(symmetry_elements.begin(),
                              symmetry_elements.end(),
                              equivalent_indices.begin(),
                              [&](const Matrix3i & g) {return ind2ind(g*r, size, centre);});

                    auto unique_indices_last = unique(equivalent_indices.begin(), equivalent_indices.end());

                    //collect measured intensities
                    auto intensities_last = transform(equivalent_indices.begin(),
                                                      unique_indices_last,
                                                      equivalent_intensities.begin(),
                                                      [&](const size_t& i){return inp.data[i];});

                    intensities_last = remove_if(equivalent_intensities.begin(), intensities_last, isnan<float>);

                    float average = average_set_of_reflections(equivalent_intensities.begin(), intensities_last, par);

                    //put the average in correct places
                    for_each(equivalent_indices.begin(),
                              unique_indices_last,
                              [&](const size_t& i){
                                  res.data[i] = average;
                                  is_reconstructed.data[i] = true;
                              });

            }


//    for(int i = 0; i < no_files; i += chunk_dims[2]) {
//        int tile_i;
//        for(tile_i = i; tile_i < min(i + static_cast<int>(chunk_dims[2]), static_cast<int>(no_files)); ++tile_i) {
//            CBFDataReader frame(argv[2 + tile_i]);
//            frame.read_data(buffer.data() + size[1]*size[2]*(tile_i - i));
//        }
//
//        hsize_t noFramesInTile = tile_i - i;
//        hsize_t memspaceOffset[] = {0, 0, 0};
//        hsize_t spaceCount[] = {noFramesInTile, memspaceSize[1], memspaceSize[2]};
//        memspace.selectHyperslab(H5S_SELECT_SET, spaceCount, memspaceOffset);
//
//        hsize_t fspaceOffset[] = {static_cast<hsize_t>(i), 0, 0};
//        fspace.selectHyperslab(H5S_SELECT_SET, spaceCount, fspaceOffset);
//
//        dataset.write(buffer.data(), PredType::NATIVE_INT32, memspace, fspace);
//    }


}

