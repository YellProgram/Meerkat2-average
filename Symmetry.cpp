//
// Created by Arkadiy on 01/06/2017.
//

#include <iostream>
#include "Symmetry.h"
#include "InputParameters.h"
#include <set>

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

    //The whole idea of generators is such that we do not need the following:
    //moreover, unique only removes duplicates which are adjacent which is useless
//    auto last = unique(res.begin(), res.end());
//    res.erase(last, res.end());

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


    if (symmetry=="m-3m") {
        return expand_generators({el(-x, y, z),
                                  el( x,-y, z),
                                  el( x, y,-z),
                                  el( y, z, x),
                                  el( y, x, z),
                                  el(-y, x, z)});
    }else if (symmetry=="m-3") {
        return expand_generators({el(-x,-y, z),
                                  el( -x,y,-z),
                                  el( z, x, y),
                                  el(-x,-y,-z)});
    }else if (symmetry=="4/mmm") {
        return expand_generators({el(-x,-y, z),
                                  el(-y, x, z),
                                  el(-x, y,-z),
                                  el(-x,-y,-z)});
    }else if (symmetry=="4/m") {
        return expand_generators({el(-x,-y, z),
                                  el(-y, x, z),
                                  el(-x,-y,-z)});
    }else if(symmetry == "mmm") {
        return expand_generators({el(-x, y, z),
                                  el( x,-y, z),
                                  el( x, y,-z)});
    }else if(symmetry == "6/mmm") {
        return expand_generators({el(-y,x-y,z),
                                  el(-x,-y,z),
                                  el(y,x,-z),
                                  el(-x,-y,-z)});
    }else if(symmetry == "6/m") {
        return expand_generators({el(-y,x-y,z),
                                  el(-x,-y,z),
                                  el(-x,-y,-z)});
    }else if(symmetry == "-3:H") {
        return expand_generators({el(-y,x-y,z),
                                 el(-x,-y,-z)});
    }else if(symmetry == "-3:R") {
        return expand_generators({el(z,x,y),
                                  el(-x,-y,-z)});
    }else if(symmetry == "-3m:H") {
        return expand_generators({el(-y,x-y,z),
                                  el( y, x,-z),
                                  el(-x,-y,-z)});
    }else if(symmetry == "-3m:R") {
        return expand_generators({el(z,x,y),
                                  el(-z,-y,-x),
                                  el(-x,-y,-z)});
    }else if(symmetry == "2/m:b") {
        return expand_generators({el(-x, y,-z),
                                  el(-x,-y,-z)});
    }else if(symmetry == "2/m:c") {
        return expand_generators({el(-x,-y, z),
                                  el(-x,-y,-z)});
    }else if(symmetry == "-1") {
        return expand_generators({el(-x, -y, -z)});
    }else if(symmetry == "1") {
        return expand_generators({el(x, y, z)});
    }else {
        assert(false);
    }
    return res;
}

///Returns linear index from the
inline size_t ind2ind(const Vector3i& r, const vector<size_t>& size,const vector<int>& centre) {
    for(int i = 0; i < 3; ++i)
        if(r[i]+centre[i] < 0 || r[i]+centre[i] >= size[i])
            return std::numeric_limits<size_t>::max();

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

inline float sum_abs(const vector<float>::iterator &start, const vector<float>::iterator &end) {
    float res = 0;
    for_each(start,
             end,
             [&](const float& i){
                 res += abs(i);
             });

    return res;
}

inline float linear_deviation_sum(const vector<float>::iterator &start, const vector<float>::iterator &end, float average) {
    float res = 0;
    for_each(start,
             end,
             [&](const float& i){
                 res += abs(i - average);
             });
    return res;
}

inline float var(const vector<float>::iterator &start, const vector<float>::iterator &end) {
    float res = 0;
    for_each(start,
             end,
             [&](const float& inp){
                 res += inp*inp;
             });
    res /= distance(start, end);
    return res;
}

inline float median(const vector<float>::iterator &start, const vector<float>::iterator &end) {
    sort(start, end);
    auto sz = distance(start, end);
    return *(start+sz/2);
}

inline const vector<float>::iterator reject_outliers(
        const vector<float>::iterator & start,
        const vector<float>::iterator & end,
        float threshold)
{
    auto sz = distance(start, end);
    float med = median(start, end);

    vector<float> differences(sz);
    transform(start, end, differences.begin(), [&](const float& i) {return abs(i - med);} );

    float median_distance = median(differences.begin(), differences.end());

    return remove_if(start, end, [&](const float& i){return abs(i - med) > median_distance*threshold;});
}



/// Returns Rint
float average(IntensityData<float>& inp, IntensityData<float>& res, const InputParameters& par) {
    vector<int> centre(3);
    for(int i=0; i<3; ++i) {
        centre[i] = round(-inp.lower_limits[i]/inp.step_sizes[i]);
    }

    auto is_reconstructed = IntensityData<bool>::empty(inp);

    auto size = inp.size;

    auto size_1d = inp.size[0]*inp.size[1]*inp.size[2];

    //Assuming the additional datasets do not exist, create them
    if(par.report_pixel_multiplicity) {
        res.other_datasets["pixel_multiplicity"] = vector<float>(size_1d, NAN);
    }
    if(par.report_pixel_rint) {
        res.other_datasets["pixel_rint"] = vector<float>(size_1d, NAN);
    }
    if(par.report_pixel_variance) {
        res.other_datasets["pixel_variance"] = vector<float>(size_1d, NAN);
    }

    double accumulated_abs_I = 0;
    double accumulated_abs_dI = 0;

    //TODO: check grid is ok with symmetry, or otherwise ignore the portion which is outside when reconstructing
    //TODO: rewrite this loop in a linear fashion??
    vector<Matrix3i> symmetry_elements = expand_symmetry(par.symmetry);
    vector<size_t> equivalent_indices(symmetry_elements.size());
    vector<float> equivalent_intensities(symmetry_elements.size());
    set<size_t> unique_indices;

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

                    //Remove elements that are outside the array. Possible in hexagonal systems and non-centered arrays
                    auto measured_ind_end = copy_if(equivalent_indices.begin(),
                                                    equivalent_indices.end(),
                                                    equivalent_indices.begin(),
                                                    [](size_t x){return x!=std::numeric_limits<size_t>::max();});

//                    auto unique_indices_last = unique(equivalent_indices.begin(), measured_ind_end);
                    unique_indices.clear();
                    unique_indices.insert(equivalent_indices.begin(), measured_ind_end);


                    //collect measured intensities
                    auto intensities_last = transform(unique_indices.begin(),
                                                      unique_indices.end(),
                                                      equivalent_intensities.begin(),
                                                      [&](const size_t& i){return inp.data[i];});

                    //remove nans
                    intensities_last = remove_if(equivalent_intensities.begin(), intensities_last, isnan<float>);

                    //reject outliers if asked
                    if(par.reject_outliers && distance(equivalent_intensities.begin(), intensities_last)>2) {
                        intensities_last = reject_outliers(equivalent_intensities.begin(), intensities_last, par.threshold);
                    }

                    float average = mean(equivalent_intensities.begin(), intensities_last);

                    float sum_abs_I = sum_abs(equivalent_intensities.begin(), intensities_last);
                    float sum_dI = linear_deviation_sum(equivalent_intensities.begin(), intensities_last, average);

                    accumulated_abs_I += sum_abs_I;
                    accumulated_abs_dI += sum_dI;

                    //put the average in correct places
                    for_each(unique_indices.begin(),
                             unique_indices.end(),
                              [&](const size_t& i){
                                  res.data[i] = average;
                                  is_reconstructed.data[i] = true;
                              });

                    //save other nonsence
                    if(par.report_pixel_multiplicity) {
                        float multiplicity = distance(equivalent_intensities.begin(), intensities_last);

                        for_each(unique_indices.begin(),
                                 unique_indices.end(),
                                 [&](const size_t& i){
                                     res.other_datasets["pixel_multiplicity"][i] = multiplicity;
                                 });
                    }

                    if(par.report_pixel_rint) {
                        for_each(unique_indices.begin(),
                                 unique_indices.end(),
                                 [&](const size_t& i){
                                     res.other_datasets["pixel_rint"][i] = sum_dI/sum_abs_I;
                                 });
                    }
                    if(par.report_pixel_variance) {
                        float variance = var(equivalent_intensities.begin(), intensities_last);

                        for_each(unique_indices.begin(),
                                 unique_indices.end(),
                                 [&](const size_t& i){
                                     res.other_datasets["pixel_variance"][i] = variance;
                                 });
                    }
                }

    return accumulated_abs_dI/accumulated_abs_I;
}

