//
// Created by Arkadiy on 01/06/2017.
//

#include <iostream>
#include "Symmetry.h"

vector<Matrix3i> expand_generators(const vector<Matrix3i>& generators) {

    auto I = Matrix3i::Identity();
    vector<Matrix3i> res = {I};

    for(auto g : generators) {
        auto last_iteration = res;
        auto gn = g;

        while(gn != I) {
            for(auto h : last_iteration)
                res.push_back(gn*h);

            gn = gn*g;
        }
    }

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
                                  el( x, y,-z)});
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

// plane_ind

void average(IntensityData<float>& inp, const string& symmetry) {

    // centre = {}
//    auto has_been_reconstructed = IntensityData<bool> inp;
    // res = new array of the same size
    //
    // for x,y,z
}