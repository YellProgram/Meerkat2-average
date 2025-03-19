//
// Created by Arkadiy on 01/06/2017.
//

#ifndef MEERKAT2_AVERAGE_SYMMETRY_H
#define MEERKAT2_AVERAGE_SYMMETRY_H
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "IntensityData.h"
#include "InputParameters.h"

using Eigen::Matrix3i;
using Eigen::Vector3i;
using namespace std;

struct RFactors {
    double R1;
    double R2;
};


vector<Matrix3i> expand_symmetry(const string& symmetry);

RFactors average(IntensityData<float>& inp, IntensityData<float>& res, const InputParameters& par);



#endif //MEERKAT2_AVERAGE_SYMMETRY_H
