#include <gtest/gtest.h>
#include <iostream>
#include "Box.hpp"

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
  EXPECT_EQ(unitbox.nDegenerate(), 0);

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

  // Get box sides
  std::map< double, std::vector< unsigned int> > sides;
  unitbox.getSides(sides);
  ASSERT_EQ(sides.size(), 1);
  double length = (sides.begin())->first;
  std::vector<unsigned int> order = (sides.begin())->second;
  EXPECT_EQ(length, 1.0);
  ASSERT_EQ(order.size(), 3);
  EXPECT_EQ(order.at(0), 0);
  EXPECT_EQ(order.at(1), 1);
  EXPECT_EQ(order.at(2), 2);

  // Get "ordered basis"
  // -> will just be order we found in
  unitbox.getBasisOrder(order);
  ASSERT_EQ(order.size(), 3);
  EXPECT_EQ(order.at(0), 0);
  EXPECT_EQ(order.at(1), 1);
  EXPECT_EQ(order.at(2), 2);

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
  EXPECT_EQ(box.nDegenerate(), 0);

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

  // Get box sides
  std::map< double, std::vector< unsigned int> > sides;
  box.getSides(sides);
  EXPECT_EQ(sides.size(), 3);
  EXPECT_TRUE(sides.find(3.0) != sides.end());
  EXPECT_TRUE(sides.find(2.0) != sides.end());
  EXPECT_TRUE(sides.find(1.0) != sides.end());

  // Get basis indices ordered in decreasing side length
  std::vector<unsigned int> order;
  box.getBasisOrder(order);
  ASSERT_EQ(order.size(), 3);
  EXPECT_EQ(order.at(0), 2);
  EXPECT_EQ(order.at(1), 1);
  EXPECT_EQ(order.at(2), 0);

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

//Tests for an oriented box
TEST(BoxTest, OrientedBox) {

  DAGMC::Vector max = {0, 1., 1.};;
  DAGMC::Vector min = {0, -1., -1.};
  DAGMC::Matrix basis = { {sqrt(2.) / 2., -sqrt(2.) / 2., 0},
    {sqrt(2.) / 2., sqrt(2.) / 2., 0},
    {0., 0., 1.}
  };


  DAGMC::Box box(max, min, basis);

  // Check construction succeeded
  ASSERT_EQ(box.getBoxStatus(), 0);
  ASSERT_TRUE(box.isSane());
  EXPECT_EQ(box.nDegenerate(), 0);

  // Get box sides
  std::map< double, std::vector< unsigned int> > sides;
  box.getSides(sides);
  EXPECT_EQ(sides.size(), 2);
  EXPECT_TRUE(sides.find(sqrt(2.)) != sides.end());
  EXPECT_TRUE(sides.find(2.) != sides.end());

  // Get basis indices ordered in decreasing side length
  std::vector<unsigned int> order;
  box.getBasisOrder(order);
  ASSERT_EQ(order.size(), 3);
  EXPECT_EQ(order.at(0), 2);
  EXPECT_EQ(order.at(1), 0);
  EXPECT_EQ(order.at(2), 1);

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


TEST(BoxTest, DegenerateBoxSquare) {

  DAGMC::Vector max = {1., 1., 0.};
  DAGMC::Vector min = {-1, -1., 0.};
  DAGMC::Box box(max, min);

  //Check construction succeeded
  ASSERT_EQ(box.getBoxStatus(), DAGMC::Box::success);
  ASSERT_TRUE(box.isSane());

  // Check degeneracy
  EXPECT_EQ(box.nDegenerate(), 1);
  EXPECT_EQ(box.getDegenDir(0), 2);

  // Get basis indices ordered in decreasing side length
  std::vector<unsigned int> order;
  box.getBasisOrder(order);
  ASSERT_EQ(order.size(), 2);
  EXPECT_EQ(order.at(0), 0);
  EXPECT_EQ(order.at(1), 1);

  // Check points
  DAGMC::Vector inside  = {0., 0., 0.};
  DAGMC::Vector outside = {0., 0., 1.};

  EXPECT_TRUE(box.containsPoint(inside));
  EXPECT_FALSE(box.containsPoint(outside));

  // Check some rays
  DAGMC::Vector orig = {0., 0., -1.};
  DAGMC::Vector dir = { 0., 0., 1,};
  EXPECT_TRUE(box.intersectsRay(orig, dir));

  dir = {1., 1., 0.};
  EXPECT_FALSE(box.intersectsRay(orig, dir));

  orig = {-2., -2., 0.};
  EXPECT_TRUE(box.intersectsRay(orig, dir));

}

TEST(BoxTest, DegenerateBoxLine) {

  DAGMC::Vector max = {1., 0., 0.};
  DAGMC::Vector min = {-1, 0., 0.};
  DAGMC::Box box(max, min);

  //Check construction succeeded
  ASSERT_EQ(box.getBoxStatus(), DAGMC::Box::success);
  ASSERT_TRUE(box.isSane());

  // Check degeneracy
  EXPECT_EQ(box.nDegenerate(), 2);
  EXPECT_EQ(box.getDegenDir(0), 1);
  EXPECT_EQ(box.getDegenDir(1), 2);

  // Get basis indices ordered in decreasing side length
  std::vector<unsigned int> order;
  box.getBasisOrder(order);
  ASSERT_EQ(order.size(), 1);
  EXPECT_EQ(order.at(0), 0);

  // Check points
  DAGMC::Vector inside  = {0., 0., 0.};
  DAGMC::Vector outside1 = {0., 0., 1.};
  DAGMC::Vector outside2 = {0., 1., 0.};
  DAGMC::Vector outside3 = {2., 0., 0.};
  EXPECT_TRUE(box.containsPoint(inside));
  EXPECT_FALSE(box.containsPoint(outside1));
  EXPECT_FALSE(box.containsPoint(outside2));
  EXPECT_FALSE(box.containsPoint(outside3));

  // Check rays

  // Intersect along line
  DAGMC::Vector orig = {-2., 0., 0.};
  DAGMC::Vector dir = { 1., 0., 0.};
  EXPECT_TRUE(box.intersectsRay(orig, dir));

  // Miss
  dir = {0, 1., 0.};
  EXPECT_FALSE(box.intersectsRay(orig, dir));
  dir = {0, 0., 1.};
  EXPECT_FALSE(box.intersectsRay(orig, dir));

  // Intersect orthogonal
  orig = {0., -1., 0.};
  dir = {0, 1., 0.};
  EXPECT_TRUE(box.intersectsRay(orig, dir));

  // Intersect at angle
  orig = {-1., -1., 0.};
  dir = {1, 1., 0.};
  EXPECT_TRUE(box.intersectsRay(orig, dir));
}

TEST(BoxTest, DegenerateBoxPoint) {

  DAGMC::Vector point = {1., 1., 1.};
  DAGMC::Box box(point, point);

  //Check construction succeeded
  ASSERT_EQ(box.getBoxStatus(), DAGMC::Box::success);
  ASSERT_TRUE(box.isSane());

  // Check degeneracy
  EXPECT_EQ(box.nDegenerate(), 3);
  EXPECT_EQ(box.getDegenDir(0), 0);
  EXPECT_EQ(box.getDegenDir(1), 1);
  EXPECT_EQ(box.getDegenDir(2), 2);

  // Get basis indices ordered in decreasing side length
  std::vector<unsigned int> order;
  box.getBasisOrder(order);
  EXPECT_EQ(order.size(), 0);

  // Check the point is contained
  EXPECT_TRUE(box.containsPoint(point));

  // Literally any other point is outside
  DAGMC::Vector outside = {0., 0., 0.};
  EXPECT_FALSE(box.containsPoint(outside));

  // Any ray passing through point should pass
  DAGMC::Vector orig = outside;
  DAGMC::Vector dir = { 1., 1., 1.};
  EXPECT_TRUE(box.intersectsRay(orig, dir));

  dir = { 1., 0., 0.};
  EXPECT_FALSE(box.intersectsRay(orig, dir));

  orig = {0., 1., 1.};
  EXPECT_TRUE(box.intersectsRay(orig, dir));

}
