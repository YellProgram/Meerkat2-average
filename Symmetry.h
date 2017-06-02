//
// Created by Arkadiy on 01/06/2017.
//

#ifndef MEERKAT2_AVERAGE_SYMMETRY_H
#define MEERKAT2_AVERAGE_SYMMETRY_H
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "IntensityData.h"

using Eigen::Matrix3i;
using Eigen::Vector3i;
using namespace std;

vector<Matrix3i> expand_symmetry(const string& symmetry);

void average(IntensityData<float>&, const string& symmetry);

#endif //MEERKAT2_AVERAGE_SYMMETRY_H
