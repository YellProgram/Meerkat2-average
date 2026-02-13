//
// Created by Arkadiy on 11/03/2016.
//

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "InputParameters.h"
#include "Exceptions.h"
#include "Symmetry.h"
#include "IntensityData.h"

using Eigen::MatrixXd;
using Eigen::Matrix;

const float eps = 0.0000001;

TEST(MeerkatAverage, Something) {
    EXPECT_EQ(1, 1);
    EXPECT_NEAR(0., 0., eps);
}

TEST(MeerkatAverage, EigenWorks) {
    Eigen::Matrix3i m;
    Eigen::Vector3i v;
    //dot product is m*m
    m(0,0) = 3;
    m(1,0) = 2;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    m(2,2) = 1;

    v << 1, 2, 3;
    m.row(0) = v;
    m << v, v, v;

    EXPECT_EQ(1, m(0,0));
}

TEST(MeerkatAverage, StringSplit) {
    auto t = split("a   b c", ' ');
    EXPECT_EQ("b", t[1]);
}

TEST(MeerkatAverage, InputParametersParser) {
    InputParameters par;
    try {
        par = parse_input("inp.mrk");
    } catch(const ParserError& e) {
        cout << e.description;
        FAIL();
    }

    EXPECT_EQ(false, par.reject_outliers);
    EXPECT_NEAR(5.7, par.threshold, eps);
    EXPECT_EQ("something_else.h5", par.input_files[1]);
    EXPECT_EQ("m-3m", par.symmetry);
    EXPECT_EQ("output.h5", par.output_name);
}

TEST(MeerkatAverage, ExpandSymmetry) {
    EXPECT_EQ(48, expand_symmetry("m-3m").size());
    EXPECT_EQ(8, expand_symmetry("mmm").size());
    EXPECT_EQ(2, expand_symmetry("-1").size());
    EXPECT_EQ(1, expand_symmetry("1").size());
    //TODO: test all other symmetry elements. At least the number of symmetry elements should be correct.
}

TEST(MeerkatAverage, DataIO) {
    auto inp = IntensityData<float>::read("test.h5", Slice());
    EXPECT_EQ(10, inp.size[0]);
    EXPECT_EQ(11, inp.size[1]);
    EXPECT_EQ(12, inp.size[2]);

    EXPECT_EQ(3, inp.lower_limits.size());

    EXPECT_NEAR(-1, inp.lower_limits.at(0), eps);
    EXPECT_NEAR(-2, inp.lower_limits.at(1), eps);
    EXPECT_NEAR(-3, inp.lower_limits.at(2), eps);

    inp.save("test-out.h5");

    inp = IntensityData<float>::read("test-out.h5", Slice());
    EXPECT_EQ(10, inp.size[0]);
    EXPECT_EQ(11, inp.size[1]);
    EXPECT_EQ(12, inp.size[2]);

    EXPECT_NEAR(-1, inp.lower_limits[0], eps);
    EXPECT_NEAR(-2, inp.lower_limits[1], eps);
    EXPECT_NEAR(-3, inp.lower_limits[2], eps);
}

TEST(MeerkatAverage, RunAveraging) {
    auto inp = IntensityData<float>::read("visual-test.h5", Slice());
    auto res = IntensityData<float>::empty(inp);
    auto par = InputParameters();
    par.symmetry = "m-3m";
    average(inp, res, InputParameters());
    res.save("visual-test-out.h5");
}
