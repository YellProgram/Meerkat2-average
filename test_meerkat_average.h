//
// Created by Arkadiy on 11/03/2016.
//

#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

#include "InputParameters.h"
#include "Exceptions.h"
#include "Symmetry.h"
#include "IntensityData.h"

using Eigen::MatrixXd;
using Eigen::Matrix;

const float eps = 0.0000001;

class MeerkatAverageTestSuit : public CxxTest::TestSuite
{
public:
    void test_something() {
        TS_ASSERT_EQUALS(1, 1);
        TS_ASSERT_DELTA(0., 0., eps);
    }

    void test_eigen_works() {
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

        TS_ASSERT_EQUALS(1, m(0,0));
    }

    void test_string_split() {
        auto t = split("a   b c", ' ');
        TS_ASSERT_EQUALS("b", t[1]);
    }

    void test_input_parameters_parser() {
        InputParameters par;
        try {
            par = parse_input("inp.mrk");
        } catch(const ParserError& e) {
            cout << e.description;
            TS_ASSERT(false);
        }

        TS_ASSERT_EQUALS(false, par.reject_outliers);
        TS_ASSERT_DELTA(5.7, par.threshold, eps);
        TS_ASSERT_EQUALS("something_else.h5", par.input_files[1]);
        TS_ASSERT_EQUALS("m3m", par.symmetry);
        TS_ASSERT_EQUALS("output.h5", par.output_name);
    }


    void test_expand_symmetry() {
        TS_ASSERT_EQUALS(48, expand_symmetry("m3m").size());
        TS_ASSERT_EQUALS(8, expand_symmetry("mmm").size());
        TS_ASSERT_EQUALS(2, expand_symmetry("-1").size());
        TS_ASSERT_EQUALS(1, expand_symmetry("1").size());
    }

    void test_data_io() {
        auto inp = IntensityData<float>::read("test.h5");
        TS_ASSERT_EQUALS(10, inp.size[0]);
        TS_ASSERT_EQUALS(11, inp.size[1]);
        TS_ASSERT_EQUALS(12, inp.size[2]);

        TS_ASSERT_DELTA(-1, inp.lower_limits[0], eps);
        TS_ASSERT_DELTA(-2, inp.lower_limits[1], eps);
        TS_ASSERT_DELTA(-3, inp.lower_limits[2], eps);

        inp.save("test-out.h5");

        inp = IntensityData<float>::read("test-out.h5");
        TS_ASSERT_EQUALS(10, inp.size[0]);
        TS_ASSERT_EQUALS(11, inp.size[1]);
        TS_ASSERT_EQUALS(12, inp.size[2]);

        TS_ASSERT_DELTA(-1, inp.lower_limits[0], eps);
        TS_ASSERT_DELTA(-2, inp.lower_limits[1], eps);
        TS_ASSERT_DELTA(-3, inp.lower_limits[2], eps);
    }
// GOT:
// matrix library eigen up and running, found matrix multiplication and vectors and matrices
// Load the instruction file
// open the instruction file, read input - source, destination, possibly several destinations? along with their scale coefficients?, symmetry, threshold, flag to kick outliers
// turn symmetry into a set of symmetry matrices
// EASY: can do this internally, by using << of eigen and only do it for mmm and m3m - the groups I care about ATM
// load reconstructions
// write out reconstructions

// TODO:

// average data without outlier rejection first

// find centre of reciprocal space index
// check that grid is compatible with the symmetry. Ignore for now hexagonal case
// check that grids of different datasets are identical


// make an array with flags to identify which data points were reconstructed
// make an array with flags on which points are outliers
// for each hkl find its equivalents
// for those select the ones which are not nans
// average or do the blessing magic: if there are more than three of those reflections, median of difference times threshold to check what is an outler, then average the ones which are not outliers
// export the averaged dataset
// possibly export the pixels which were averaged

};

