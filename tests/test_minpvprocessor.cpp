/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <config.h>

#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE MinpvProcessorTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/MinpvProcessor.hpp>


BOOST_AUTO_TEST_CASE(GAP_MAXGAP)
{
    // Set up a simple example.
    std::vector<double> zcorn = { 0, 0, 0, 0,
                                  2, 2, 2, 2,
                                  2, 2, 2, 2,
                                  2.5, 2.5, 2.5, 2.5,
                                  2.8, 2.8, 2.8, 2.8,
                                  3.5, 3.5, 3.5, 3.5};
    std::vector<double> pv = { 2, 0.5, 0.7};
    std::vector<double> minpvv(3, 0.6);
    std::vector<int> actnum = { 1, 1, 1 };
    std::vector<double> thickness = {2, 0.5, 0.7};
    double z_threshold = 0.4;

    Opm::MinpvProcessor mp1(1, 1, 3);
    auto z1 = zcorn;
    double max_gap = 1e20;
    bool fill_removed_cells = false;
    bool pinch_no_gap = false;

    auto minpv_result = mp1.process(thickness, z_threshold, max_gap, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);
    BOOST_CHECK_EQUAL(minpv_result.nnc[0], 2);

    max_gap = .29;
    minpv_result = mp1.process(thickness, z_threshold, max_gap, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 0);
}

BOOST_AUTO_TEST_CASE(GAP_MAXGAP_no_pinched_cells)
{
    // Set up a simple example.
    std::vector<double> zcorn = { 0, 0, 0, 0,
                                  2, 2.1, 2, 2,
                                  2, 2.1, 2, 2,
                                  2.5, 2.5, 2.5, 2.5,
                                  2.8, 2.8, 2.8, 2.8,
                                  3.5, 3.5, 3.5, 3.5};
    std::vector<double> pv = { 2, 0.5, 0.7};
    std::vector<double> minpvv(3, 0.0);
    std::vector<int> actnum = { 1, 1, 1 };
    std::vector<double> thickness = {2, 0.5, 0.7};
    double z_threshold = 0.0;

    Opm::MinpvProcessor mp1(1, 1, 3);
    double max_gap = 1e20;
    bool fill_removed_cells = false;
    bool pinch_no_gap = false;

    // Use options that will create NNCs for vertically unconnected cells with small gaps without cells being pinched.
    auto minpv_result = mp1.process(thickness, z_threshold, max_gap, pv, minpvv, actnum, fill_removed_cells,
                                    zcorn.data(), pinch_no_gap, false, {}, [](int){ return 1; });
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);
    if (minpv_result.nnc.size() )
      BOOST_CHECK_EQUAL(minpv_result.nnc[1], 2);
}

BOOST_AUTO_TEST_CASE(Pinch4ALL)
{
    // Set up a simple example.
    std::vector<double> zcorn = { 0, 0, 0, 0,
                                  2, 2, 2, 2,
                                  2, 2, 2, 2,
                                  2.5, 2.5, 2.5, 2.5,
                                  2.5, 2.5, 2.5, 2.5,
                                  3.5, 3.5, 3.5, 3.5,
                                  3.5, 3.5, 3.5, 3.5,
                                  6, 6, 6, 6 };

    std::vector<double> pv = { 2, 0.5, 0.5, 2.5};
    std::vector<int> actnum = { 1, 1, 1, 1 };
    std::vector<double> thickness = {2, 0.5, 0.5, 2.5};
    std::vector<double> permz = {2, 2, 0, 2.5};
    auto multz = [](int){ return 1.0;};
    double z_threshold = 0.4;

    Opm::MinpvProcessor mp1(1, 1, 4);
    auto z1 = zcorn;
    std::vector<double> minpvv(4, 0.6);
    bool fill_removed_cells = false;
    bool pinch_no_gap = false;
    bool option4all = true;

    // Test PINCH option 4 being ALL
    auto minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum,
				    fill_removed_cells, z1.data(), pinch_no_gap,
				    option4all, permz, multz);
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 0);

    permz = {2, 2, 1, 2.5}; // Should create an NNC
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum,
                               fill_removed_cells, z1.data(), pinch_no_gap,
                               option4all, permz, multz);
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);
    auto multz2 = [](int i){ if (i==2) return 0.0; else return 1.0;};

    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum,
			       fill_removed_cells, z1.data(), pinch_no_gap,
			       option4all, permz, multz2);
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 0);

    auto multz3 = [](int i){ if (i==0) return 0.0; else return 1.0;};
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum,
			       fill_removed_cells, z1.data(), pinch_no_gap,
			       option4all, permz, multz3);
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 0);
}
BOOST_AUTO_TEST_CASE(Pinch)
{
    // Set up a simple example.
    std::vector<double> zcorn = { 0, 0, 0, 0,
                                  2, 2, 2, 2,
                                  2, 2, 2, 2,
                                  2.5, 2.5, 2.5, 2.5,
                                  2.5, 2.5, 2.5, 2.5,
                                  3.5, 3.5, 3.5, 3.5,
                                  3.5, 3.5, 3.5, 3.5,
                                  6, 6, 6, 6 };

    std::vector<double> pv = { 2, 0.5, 1, 2.5};
    std::vector<int> actnum = { 1, 1, 1, 1 };
    std::vector<double> thickness = {2, 0.5, 1, 2.5};
    double z_threshold = 0.4;

    Opm::MinpvProcessor mp1(1, 1, 4);
    auto z1 = zcorn;
    std::vector<double> minpvv(4, 0.6);
    bool fill_removed_cells = false;
    bool pinch_no_gap = false;
    std::vector<double> zcornAfter =
        {0, 0, 0, 0,
         2, 2, 2, 2,
         2, 2, 2, 2,
         2, 2, 2, 2,
         2.5, 2.5, 2.5, 2.5,
         3.5, 3.5, 3.5, 3.5,
         3.5, 3.5, 3.5, 3.5,
         6, 6, 6, 6
        };
    auto minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcornAfter.begin(), zcornAfter.end());
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);

    z1= zcorn;
    pinch_no_gap = true;
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcornAfter.begin(), zcornAfter.end());
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 0); // No nnc because deactivated cell thickness is too high

    z_threshold = 0.51;
    pinch_no_gap = true;
    minpvv = std::vector(4, 0.6);
    z1 = zcorn;
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);

    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);
    BOOST_CHECK_EQUAL(minpv_result.nnc[0], 2);
    BOOST_CHECK(minpv_result.removed_cells == std::vector<std::size_t>{1});
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcornAfter.begin(), zcornAfter.end());

    z_threshold = 1.1;
    pinch_no_gap = true;
    minpvv = std::vector(4, 1.1);
    z1 = zcorn;
    double max_gap = 1.0;
    minpv_result = mp1.process(thickness, z_threshold, max_gap, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);

    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 0);
    BOOST_CHECK((minpv_result.removed_cells == std::vector<std::size_t>{1, 2}));

    z_threshold = 0.4;
    pinch_no_gap = false;
    minpvv = std::vector(4, 0.6);
    z1 = zcorn;
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);

    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);
    BOOST_CHECK(minpv_result.removed_cells == std::vector<std::size_t>{1});
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcornAfter.begin(), zcornAfter.end());

    z_threshold = 0.6;
    pinch_no_gap = true;
    minpvv = std::vector(4, 0.4);
    z1 = zcorn;
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);

    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcorn.begin(), zcorn.end());

    z_threshold = 1.1;
    pinch_no_gap = false;
    minpvv = std::vector(4, 1.1);
    z1 = zcorn;
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);
    zcornAfter =
        {0, 0, 0, 0,
         2, 2, 2, 2,
         2, 2, 2, 2,
         2, 2, 2, 2,
         2.5, 2.5, 2.5, 2.5,
         2.5, 2.5, 2.5, 2.5,
         3.5, 3.5, 3.5, 3.5,
         6, 6, 6, 6
        };
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcornAfter.begin(), zcornAfter.end());

    z_threshold = 1.1;
    pinch_no_gap = true;
    minpvv = std::vector(4, 1.1);
    z1 = zcorn;
    minpv_result = mp1.process(thickness, z_threshold, 1e20, pv, minpvv, actnum, fill_removed_cells, z1.data(), pinch_no_gap);
    zcornAfter =
        {0, 0, 0, 0,
         2, 2, 2, 2,
         2, 2, 2, 2,
         2, 2, 2, 2,
         2.5, 2.5, 2.5, 2.5,
         2.5, 2.5, 2.5, 2.5,
         3.5, 3.5, 3.5, 3.5,
         6, 6, 6, 6
        };
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1);
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcornAfter.begin(), zcornAfter.end());
    auto exp =  std::vector<std::size_t>{1, 2};
    BOOST_CHECK(minpv_result.removed_cells == exp);
    BOOST_CHECK_EQUAL(minpv_result.nnc.size(), 1u);
    BOOST_CHECK_EQUAL(minpv_result.nnc[0], 3);
}

BOOST_AUTO_TEST_CASE(Processing)
{
    std::vector<double> zcorn = { 0, 0, 0, 0,
                                  2, 2, 2, 2,
                                  2, 2, 2, 2,
                                  3, 3, 3, 3,
                                  3, 3, 3, 3,
                                  3.1, 3.1, 3.1, 3.1,
                                  3.1, 3.1, 3.1, 3.1,
                                  6, 6, 6, 6 };
    std::vector<double> zcorn1after = { 0, 0, 0, 0,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        3, 3, 3, 3,
                                        3, 3, 3, 3,
                                        3, 3, 3, 3, // thin inactive cell collapsed
                                        3.1, 3.1, 3.1, 3.1,
                                        6, 6, 6, 6 };
    std::vector<double> zcorn1bafter = { 0, 0, 0, 0,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        3, 3, 3, 3,
                                        3, 3, 3, 3,
                                        3, 3, 3, 3, // collapsed thin active cell with small pv
                                        3, 3, 3, 3, // pv filled,
                                        6, 6, 6, 6 };
    std::vector<double> zcorn2after = { 0, 0, 0, 0,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        6, 6, 6, 6 };
    std::vector<double> zcorn3after = { 0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        6, 6, 6, 6  };
    std::vector<double> zcorn4after = { 0, 0, 0, 0,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        3.1, 3.1, 3.1, 3.1,
                                        6, 6, 6, 6 };
    std::vector<double> zcorn5after = { 0, 0, 0, 0,
                                        0, 0, 0, 0,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        3.1, 3.1, 3.1, 3.1,
                                        6, 6, 6, 6 };
    std::vector<double> zcorn7after = { 0, 0, 0, 0,
                                        2, 2, 2, 2,
                                        2, 2, 2, 2,
                                        3, 3, 3, 3,
                                        3, 3, 3, 3,
                                        3, 3, 3, 3,
                                        3.1, 3.1, 3.1, 3.1,
                                        6, 6, 6, 6 };

    std::vector<double> pv = { 2, 1, 0, 3};
    std::vector<int> actnum = { 1, 1, 0, 1 };
    std::vector<int> actnum_empty;
    std::vector<double> thicknes = {2, 1, 0, 3};
    double z_threshold = 0.0;

    Opm::MinpvProcessor mp1(1, 1, 4);
    auto z1 = zcorn;
    std::vector<double> minpvv1(4, 0.5);
    bool fill_removed_cells = true;
    mp1.process(thicknes, z_threshold, 1e20, pv, minpvv1, actnum, fill_removed_cells, z1.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z1.begin(), z1.end(), zcorn1after.begin(), zcorn1after.end());

    // check that it is possible to pass empty actnum
    Opm::MinpvProcessor mp1b(1, 1, 4);
    auto z1b = zcorn;
    mp1b.process(thicknes, z_threshold, 1e20, pv, minpvv1, actnum_empty, fill_removed_cells, z1b.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z1b.begin(), z1b.end(), zcorn1bafter.begin(), zcorn1bafter.end());

    Opm::MinpvProcessor mp2(1, 1, 4);
    auto z2 = zcorn;
    std::vector<double> minpvv2(4, 1.5);
    mp2.process(thicknes, z_threshold, 1e20, pv, minpvv2, actnum, fill_removed_cells, z2.data());
    BOOST_CHECK_EQUAL_COLLECTIONS(z2.begin(), z2.end(), zcorn2after.begin(), zcorn2after.end());

    Opm::MinpvProcessor mp3(1, 1, 4);
    auto z3 = zcorn;
    std::vector<double> minpvv3(4, 2.5);
    auto minpv_result3 = mp3.process(thicknes, z_threshold, 1e20, pv, minpvv3, actnum, fill_removed_cells, z3.data());
    BOOST_CHECK_EQUAL(minpv_result3.nnc.size(), 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(z3.begin(), z3.end(), zcorn3after.begin(), zcorn3after.end());

    Opm::MinpvProcessor mp4(1, 1, 4);
    auto z4 = zcorn;
    auto minpv_result4 = mp4.process(thicknes, z_threshold, 1e20, pv, minpvv2, actnum, !fill_removed_cells, z4.data());
    BOOST_CHECK_EQUAL(minpv_result4.nnc.size(), 1);
    BOOST_CHECK_EQUAL(minpv_result4.nnc.at(0), 3);
    BOOST_CHECK(minpv_result4.removed_cells == std::vector<std::size_t>{1});
    BOOST_CHECK_EQUAL_COLLECTIONS(z4.begin(), z4.end(), zcorn4after.begin(), zcorn4after.end());

    Opm::MinpvProcessor mp5(1, 1, 4);
    auto z5 = zcorn;
    auto minpv_result5 = mp5.process(thicknes, z_threshold, 1e20, pv, minpvv3, actnum, !fill_removed_cells, z5.data());
    BOOST_CHECK_EQUAL(minpv_result5.nnc.size(), 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(z5.begin(), z5.end(), zcorn5after.begin(), zcorn5after.end());

    Opm::MinpvProcessor mp6(1, 1, 4);
    auto z6 = zcorn;
    std::vector<double> minpvv6 = {1, 2, 2, 1};
    auto minpv_result6 = mp6.process(thicknes, z_threshold, 1e20, pv, minpvv6, actnum, !fill_removed_cells, z6.data());
    BOOST_CHECK_EQUAL(minpv_result6.nnc.size(), 1);
    BOOST_CHECK_EQUAL_COLLECTIONS(z6.begin(), z6.end(), zcorn4after.begin(), zcorn4after.end());

    Opm::MinpvProcessor mp7(1, 1, 4);
    auto z7 = zcorn;
    std::vector<double> minpvv7(4, 0); // don't deactivate any cells
    thicknes = {2, 1, 0.1, 3};
    z_threshold = 0.2; // create NNC over it
    auto minpv_result7 = mp7.process(thicknes, z_threshold, 1e20, pv, minpvv7, actnum, !fill_removed_cells, z7.data());
    BOOST_CHECK_EQUAL(minpv_result7.nnc.size(), 1);
    BOOST_CHECK_EQUAL_COLLECTIONS(z7.begin(), z7.end(), zcorn7after.begin(), zcorn7after.end());
}
