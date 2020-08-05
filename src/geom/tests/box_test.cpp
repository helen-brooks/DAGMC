#include <gtest/gtest.h>
#include <iostream>
//#include "libmesh/libmesh.h"
#include "box.hpp"

// Box test suite
// Tests organised by box type rather than method

// Tests for a unit box
TEST(BoxTest, UnitBox) {

  // Build a unit box
  DAGMC::Vector maxpoint = {1., 1., 1.};
  DAGMC::Box unitbox(maxpoint);

  // Check construction succeeded
  ASSERT_EQ(unitbox.getBoxStatus(), 0);
  ASSERT_TRUE(unitbox.isSane());

  // Check basis vectors components
  DAGMC::Vector v1 = unitbox.basisVec(0);
  EXPECT_EQ(v1(0), 1.);
  EXPECT_EQ(v1(1), 0.);
  EXPECT_EQ(v1(2), 0.);

  DAGMC::Vector v2 = unitbox.basisVec(1);
  EXPECT_EQ(v2(0), 0.);
  EXPECT_EQ(v2(1), 1.);
  EXPECT_EQ(v2(2), 0.);

  DAGMC::Vector v3 = unitbox.basisVec(2);
  EXPECT_EQ(v3(0), 0.);
  EXPECT_EQ(v3(1), 0.);
  EXPECT_EQ(v3(2), 1.);

  // Try to fetch an out-of-bounds basis vector
  EXPECT_THROW(unitbox.basisVec(3), std::logic_error);

  // Check some points
  DAGMC::Vector inside1 = {0.5, 0.5, 0.5};
  DAGMC::Vector inside2 = {0.25, 0.75, 0.1};
  DAGMC::Vector inside3 = {0.75, 0.1, 0.25};
  DAGMC::Vector inside4 = {0.1, 0.25, 0.75};
  EXPECT_TRUE(unitbox.containsPoint(inside1));
  EXPECT_TRUE(unitbox.containsPoint(inside2));
  EXPECT_TRUE(unitbox.containsPoint(inside3));
  EXPECT_TRUE(unitbox.containsPoint(inside4));

  DAGMC::Vector outside1 = {2.0, 1.0, 1.0};
  DAGMC::Vector outside2 = {1.0, 2.0, 1.0};
  DAGMC::Vector outside3 = {1.0, 1.0, 2.0};
  DAGMC::Vector outside4 = {-1., -1., -1.};
  EXPECT_FALSE(unitbox.containsPoint(outside1));
  EXPECT_FALSE(unitbox.containsPoint(outside2));
  EXPECT_FALSE(unitbox.containsPoint(outside3));
  EXPECT_FALSE(unitbox.containsPoint(outside4));

  // Edge cases should still be in box
  DAGMC::Vector edge1 = {0.5, 0., 0.};
  DAGMC::Vector edge2 = {0., 0.5, 0.};
  DAGMC::Vector edge3 = {0., 0., 0.5};
  EXPECT_TRUE(unitbox.containsPoint(edge1));
  EXPECT_TRUE(unitbox.containsPoint(edge2));
  EXPECT_TRUE(unitbox.containsPoint(edge3));

  // Corners should be in box
  EXPECT_TRUE(unitbox.containsPoint(v1));
  EXPECT_TRUE(unitbox.containsPoint(v2));
  EXPECT_TRUE(unitbox.containsPoint(v3));
  EXPECT_TRUE(unitbox.containsPoint(maxpoint));

  // Check some ray intersections

  // Ray along x
  DAGMC::Vector orig1 = {-1.0, 0.5, 0.5};
  DAGMC::Vector orig2 = { 2.0, 0.5, 0.5};
  double tIn(0.);
  double tOut(0.);
  //Box in front
  EXPECT_TRUE(unitbox.intersectsRay(orig1, v1, tIn, tOut));
  EXPECT_GT(tIn, 0.);
  EXPECT_GT(tOut, 0.);
  EXPECT_TRUE(unitbox.intersectsRay(orig2, -v1, tIn, tOut));
  EXPECT_GT(tIn, 0.);
  EXPECT_GT(tOut, 0.);

  //Box behind
  EXPECT_TRUE(unitbox.intersectsRay(orig1, -v1, tIn, tOut));
  EXPECT_LT(tIn, 0.);
  EXPECT_LT(tOut, 0.);
  EXPECT_TRUE(unitbox.intersectsRay(orig1, v1));
  EXPECT_LT(tIn, 0.);
  EXPECT_LT(tOut, 0.);

  // Ray along y
  orig1 = {0.5, -1.0, 0.5};
  orig2 = {0.5, 2.0, 0.5};
  EXPECT_TRUE(unitbox.intersectsRay(orig1, v2));
  EXPECT_TRUE(unitbox.intersectsRay(orig2, -v2));

  // Ray along z
  orig1 = {0.5, 0.5, -1.0};
  orig2 = {0.5, 0.5,  2.0};
  EXPECT_TRUE(unitbox.intersectsRay(orig1, v3));
  EXPECT_TRUE(unitbox.intersectsRay(orig2, -v3));

  // Ray at angle
  DAGMC::Vector orig3 = {-1., 0, 0};
  DAGMC:: Vector dir1 = {1., 0.5, 0.};
  DAGMC:: Vector dir2 = {1., 2.0, 0.};
  EXPECT_TRUE(unitbox.intersectsRay(orig3, dir1));
  EXPECT_FALSE(unitbox.intersectsRay(orig3, dir2));

  // Check some misses
  EXPECT_FALSE(unitbox.intersectsRay(orig3, v2));
  EXPECT_FALSE(unitbox.intersectsRay(orig3, v3));

  // Some edge cases:
  // Intersect along side
  EXPECT_TRUE(unitbox.intersectsRay(orig3, v1));

  // Intersect at corner
  dir1 = {1., 1., 0.};
  dir2 = {1., 0., 1.};
  EXPECT_TRUE(unitbox.intersectsRay(orig3, dir1));
  EXPECT_TRUE(unitbox.intersectsRay(orig3, dir2));

}

//Tests for an aligned box
TEST(BoxTest, AlignedBox) {

  // Build a aligned box 1 x 2 x 3
  DAGMC::Vector maxpoint = {2., 3., 4.};
  DAGMC::Vector minpoint = {1., 1., 1.};

  DAGMC::Box box(maxpoint, minpoint);

  // Check construction succeeded
  ASSERT_EQ(box.getBoxStatus(), 0);
  ASSERT_TRUE(box.isSane());

  // Check basis vectors components
  DAGMC::Vector v1 = box.basisVec(0);
  EXPECT_EQ(v1(0), 1.);
  EXPECT_EQ(v1(1), 0.);
  EXPECT_EQ(v1(2), 0.);

  DAGMC::Vector v2 = box.basisVec(1);
  EXPECT_EQ(v2(0), 0.);
  EXPECT_EQ(v2(1), 1.);
  EXPECT_EQ(v2(2), 0.);

  DAGMC::Vector v3 = box.basisVec(2);
  EXPECT_EQ(v3(0), 0.);
  EXPECT_EQ(v3(1), 0.);
  EXPECT_EQ(v3(2), 1.);

  // Check some points
  DAGMC::Vector inside1 = {1.5, 2.0, 2.0};
  DAGMC::Vector inside2 = {1.5, 2.0, 3.0};
  DAGMC::Vector inside3 = {1.25, 1.25, 1.25};
  DAGMC::Vector inside4 = {1.95, 2.95, 3.95};
  EXPECT_TRUE(box.containsPoint(inside1));
  EXPECT_TRUE(box.containsPoint(inside2));
  EXPECT_TRUE(box.containsPoint(inside3));
  EXPECT_TRUE(box.containsPoint(inside4));

  DAGMC::Vector outside1 = {-1., -1, -1};
  DAGMC::Vector outside2 = {0., 0., 0.};
  DAGMC::Vector outside3 = {2., 3., 5.};
  DAGMC::Vector outside4 = {1.5, 2.0, 5.0};
  EXPECT_FALSE(box.containsPoint(outside1));
  EXPECT_FALSE(box.containsPoint(outside2));
  EXPECT_FALSE(box.containsPoint(outside3));
  EXPECT_FALSE(box.containsPoint(outside4));

  // Check some edge points
  DAGMC::Vector edge1 = {1.5, 1., 1.};
  DAGMC::Vector edge2 = {1., 2., 1.};
  DAGMC::Vector edge3 = {1., 1., 3.};
  DAGMC::Vector edge4 = {1., 2., 3.};
  EXPECT_TRUE(box.containsPoint(edge1));
  EXPECT_TRUE(box.containsPoint(edge1));
  EXPECT_TRUE(box.containsPoint(edge1));
  EXPECT_TRUE(box.containsPoint(edge1));

  // Check some corners
  DAGMC::Vector corner1 = {2., 1., 1.};
  DAGMC::Vector corner2 = {1., 3., 1.};
  DAGMC::Vector corner3 = {1., 1., 4.};
  DAGMC::Vector corner4 = {1., 3., 4.};
  EXPECT_TRUE(box.containsPoint(minpoint));
  EXPECT_TRUE(box.containsPoint(maxpoint));
  EXPECT_TRUE(box.containsPoint(corner1));
  EXPECT_TRUE(box.containsPoint(corner2));
  EXPECT_TRUE(box.containsPoint(corner3));
  EXPECT_TRUE(box.containsPoint(corner4));

  // Check some ray intersections

  // Ray along x
  DAGMC::Vector orig1 = {0., 1., 1.};
  DAGMC::Vector orig2 = {0., 2., 3.};
  DAGMC::Vector orig3 = {0., 3., 3.};
  EXPECT_TRUE(box.intersectsRay(orig1, v1));
  EXPECT_TRUE(box.intersectsRay(orig1, -v1));
  EXPECT_TRUE(box.intersectsRay(orig2, v1));
  EXPECT_TRUE(box.intersectsRay(orig3, v1));

  // Ray along y
  EXPECT_FALSE(box.intersectsRay(orig1, v2));
  orig1 = {1., 0., 1.};
  orig2 = {1.5, 0., 3.};
  orig3 = {2., 0., 4.};
  EXPECT_TRUE(box.intersectsRay(orig1, v2));
  EXPECT_TRUE(box.intersectsRay(orig1, -v2));
  EXPECT_TRUE(box.intersectsRay(orig2, v2));
  EXPECT_TRUE(box.intersectsRay(orig3, v2));

  // Ray along z
  EXPECT_FALSE(box.intersectsRay(orig1, v3));
  orig1 = {1., 1., 0.};
  orig2 = {1.5, 2., 0.};
  orig3 = {2., 3., 0.};
  EXPECT_TRUE(box.intersectsRay(orig1, v3));
  EXPECT_TRUE(box.intersectsRay(orig1, -v3));
  EXPECT_TRUE(box.intersectsRay(orig2, v3));
  EXPECT_TRUE(box.intersectsRay(orig3, v3));

  // Ray on digaonals
  DAGMC::Vector d1 = {1., 2., 3.};
  DAGMC::Vector d2 = {-1., 2., 3.};
  orig1 = {3., 0., 0.};
  orig2 = {1.95, 1.5, 1.};
  orig3 = {1.25, 2.95, 1.};
  orig2 = orig2 - d1;
  orig3 = orig3 - d1;

  EXPECT_FALSE(box.intersectsRay(orig1, d1));
  EXPECT_TRUE(box.intersectsRay(orig2, d1));
  EXPECT_TRUE(box.intersectsRay(orig3, d1));
  EXPECT_TRUE(box.intersectsRay(corner1, d1));
  EXPECT_TRUE(box.intersectsRay(minpoint, d1));

  orig1 = {1.5, 1.5, 1.5};
  orig2 = orig1 - d1;
  orig3 = corner4 - d2;
  EXPECT_TRUE(box.intersectsRay(orig1, d2));
  EXPECT_FALSE(box.intersectsRay(orig2, d2));
  EXPECT_TRUE(box.intersectsRay(orig3, d2));
  EXPECT_TRUE(box.intersectsRay(corner4, d2));
  EXPECT_TRUE(box.intersectsRay(maxpoint, d2));

}

// Test for boxes where dimensionalities don't match
TEST(BoxTest, BrokenBoxDim) {

  DAGMC::Vector max = {1., 1., 1.};
  DAGMC::Vector min = {0., 0.};
  DAGMC::Box box1(max, min);

  EXPECT_EQ(box1.getBoxStatus(), DAGMC::Box::faildim);
  EXPECT_FALSE(box1.isSane());

  min = {0., 0., 0.};
  DAGMC::Matrix id = { {1., 0.,}, {0., 1.} };
  DAGMC::Box box2(max, min, id);
  EXPECT_EQ(box2.getBoxStatus(), DAGMC::Box::faildim);
  EXPECT_FALSE(box2.isSane());

}

// Test for boxes where the eigenvectors are not normalised
TEST(BoxTest, BrokenBoxNorm) {

  DAGMC::Vector max = {1., 1., 1.};
  DAGMC::Vector min = {0., 0., 0.};
  DAGMC::Matrix basis1 = { {2., 0., 0}, {0., 1., 0.}, {0., 0., 1.} };
  DAGMC::Matrix basis2 = { {1., 0., 0}, {0., 2., 0.}, {0., 0., 1.} };
  DAGMC::Matrix basis3 = { {1., 0., 0}, {0., 1., 0.}, {0., 0., 2.} };

  DAGMC::Box box1(max, min, basis1);
  DAGMC::Box box2(max, min, basis2);
  DAGMC::Box box3(max, min, basis3);

  EXPECT_EQ(box1.getBoxStatus(), DAGMC::Box::failnorm);
  EXPECT_FALSE(box1.isSane());
  EXPECT_EQ(box2.getBoxStatus(), DAGMC::Box::failnorm);
  EXPECT_FALSE(box2.isSane());
  EXPECT_EQ(box3.getBoxStatus(), DAGMC::Box::failnorm);
  EXPECT_FALSE(box3.isSane());

}

// Test for boxes where the eigenvectors are not orthogonal
TEST(BoxTest, BrokenBoxOrth) {

  DAGMC::Vector max = {1., 1., 1.};
  DAGMC::Vector min = {0., 0., 0.};

  DAGMC::Matrix basis1 = { {1., sqrt(2) / 2., 0}, {0., sqrt(2) / 2., 0.}, {0., 0., 1.} };
  DAGMC::Matrix basis2 = { {1., 0., 0}, {0., 1, sqrt(2) / 2.}, {0., 0., sqrt(2) / 2.} };
  DAGMC::Matrix basis3 = { {1., 0., sqrt(2) / 2.}, {0., 1, 0.}, {0., 0., sqrt(2) / 2.} };

  DAGMC::Box box1(max, min, basis1);
  DAGMC::Box box2(max, min, basis2);
  DAGMC::Box box3(max, min, basis3);

  EXPECT_EQ(box1.getBoxStatus(), DAGMC::Box::failorth);
  EXPECT_FALSE(box1.isSane());
  EXPECT_EQ(box2.getBoxStatus(), DAGMC::Box::failorth);
  EXPECT_FALSE(box2.isSane());
  EXPECT_EQ(box3.getBoxStatus(), DAGMC::Box::failorth);
  EXPECT_FALSE(box3.isSane());

}

TEST(BoxTest, BrokenBoxOrder) {

  DAGMC::Vector max = {1., 1., 1.};
  DAGMC::Vector min = {0., 0., 0.};
  DAGMC::Matrix basis1 = { {1., 0., 0}, {0., 1., 0.}, {0., 0., 1.} };
  DAGMC::Matrix basis2 = -1.*basis1;

  // Min / max wrong way round
  DAGMC::Box box1(min, max, basis1);

  // Min/ max "correct" ordering, but basis is reversed
  DAGMC::Box box2(max, min, basis2);

  EXPECT_EQ(box1.getBoxStatus(), DAGMC::Box::failordered);
  EXPECT_FALSE(box1.isSane());
  EXPECT_EQ(box2.getBoxStatus(), DAGMC::Box::failordered);
  EXPECT_FALSE(box2.isSane());

}

TEST(BoxTest, BrokenBoxDegenerate) {

  DAGMC::Vector max = {0., 1., 0.};
  DAGMC::Vector min = {0, -1., 0.};
  DAGMC::Box box(min, max);

  EXPECT_EQ(box.getBoxStatus(), DAGMC::Box::faildegenerate);
  EXPECT_FALSE(box.isSane());

}

//Tests for an oriented box
TEST(BoxTest, OrientedBox) {

  DAGMC::Vector max = {0, 1., 1.};;
  DAGMC::Vector min = {0, -1., -1.};
  DAGMC::Matrix basis = { {sqrt(2) / 2., -sqrt(2) / 2., 0},
    {sqrt(2) / 2., sqrt(2) / 2., 0},
    {0., 0., 1.}
  };


  DAGMC::Box box(max, min, basis);

  // Check construction succeeded
  ASSERT_EQ(box.getBoxStatus(), 0);
  ASSERT_TRUE(box.isSane());

  // Check some points
  DAGMC::Vector orig1 = { 0., 0., 0.};
  DAGMC::Vector orig2 = {-1., -1., 0.};
  DAGMC::Vector orig3 = { 1., 1., 0.};
  DAGMC::Vector orig4 = {-1,  0., 0.};
  DAGMC::Vector orig5 = { 1., 0., 0.};
  DAGMC::Vector orig6 = {-2., 0., 0.};

  EXPECT_TRUE(box.containsPoint(orig1));
  EXPECT_FALSE(box.containsPoint(orig2));
  EXPECT_FALSE(box.containsPoint(orig3));
  EXPECT_TRUE(box.containsPoint(orig4));
  EXPECT_TRUE(box.containsPoint(orig5));
  EXPECT_FALSE(box.containsPoint(orig6));

  // Check some rays
  DAGMC::Vector dir1 = {sqrt(2) / 2., sqrt(2) / 2., 0};
  DAGMC::Vector dir2 = {0., 1., 0.};
  DAGMC::Vector dir3 = {0., 0., 1.};

  EXPECT_TRUE(box.intersectsRay(orig1, dir1));
  EXPECT_TRUE(box.intersectsRay(orig1, dir2));
  EXPECT_TRUE(box.intersectsRay(orig1, dir3));
  EXPECT_TRUE(box.intersectsRay(orig2, dir1));
  EXPECT_TRUE(box.intersectsRay(orig2, dir2));
  EXPECT_FALSE(box.intersectsRay(orig2, dir3));
  EXPECT_TRUE(box.intersectsRay(orig3, dir1));
  EXPECT_TRUE(box.intersectsRay(orig4, dir1));
  EXPECT_TRUE(box.intersectsRay(orig5, dir1));
  EXPECT_FALSE(box.intersectsRay(orig6, dir1));

}
