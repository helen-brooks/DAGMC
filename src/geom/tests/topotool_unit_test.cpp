#include <gtest/gtest.h>
#include <iostream>
#include "libmesh/libmesh.h"
#include "GeomTopoTool.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class GeomTopoToolLMTest : public ::testing::Test {
 protected:

  GeomTopoToolLMTest() : meshLoaded(false),
    libMeshException(false) {};

  // Initalize variables for each test
  virtual void SetUp() {

    //Emulate dummy command line args since required by libmesh
    int argc_dummy = 1;
    char dummych[] = "dummy";
    char* argv_dummy[] = { dummych };

    // Initialize libMesh and handle any exceptions
    libmesh_try {
      // Initialize the library
      initPtr = std::make_shared<libMesh::LibMeshInit>(argc_dummy, argv_dummy);
      // Create the mesh
      if (initPtr != nullptr) {
        meshPtr = std::make_shared<libMesh::Mesh>(initPtr->comm());
      }

      // Set the input file
      std::string infile = "sphere.e";

      // Read in the mesh data
      if (meshPtr != nullptr) {
        meshPtr->read(infile);
        meshLoaded = meshPtr->is_prepared();
      }
    }
    libmesh_catch(libMesh::LogicError & e) {
      libMeshException = true;
      return;
    }

    // Create the geomtopotool
    gTTPtr = std::make_shared<DAGMC::GeomTopoToolLM>(meshPtr);

  }

  // Nothing to do
  virtual void TearDown() {}

  // Data members required for tests
  std::shared_ptr<libMesh::LibMeshInit> initPtr;
  std::shared_ptr<libMesh::Mesh> meshPtr;
  std::shared_ptr<DAGMC::GeomTopoTool> gTTPtr;
  bool meshLoaded;
  bool libMeshException;

};

//---------------------------------------------------------------------------//
// FIXTURE BASED TESTS
//---------------------------------------------------------------------------//

// Test that we set  up correctly
TEST_F(GeomTopoToolLMTest, SetUp) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);

  EXPECT_NE(initPtr, nullptr);
  EXPECT_NE(meshPtr, nullptr);
  EXPECT_NE(gTTPtr, nullptr);
  EXPECT_TRUE(meshLoaded);

}

// Test that we can set up geometry
TEST_F(GeomTopoToolLMTest, setupGeometry) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);

  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_NE(gTTPtr, nullptr);
  ASSERT_TRUE(meshLoaded);

  // Check we don't fail setup
  ASSERT_TRUE(gTTPtr->setupGeometry());

  // Check we get the right numbers of surfs/vols
  EXPECT_EQ(gTTPtr->nVols(), 2);
  EXPECT_EQ(gTTPtr->nSurfs(), 1);

  // Check volume indices
  EXPECT_EQ(gTTPtr->getVolID(0), DAGMC::IMPLICIT_COMPLEMENT);
  EXPECT_EQ(gTTPtr->getVolID(1), 1);
  EXPECT_EQ(gTTPtr->getVolID(2), DAGMC::VOID_INDEX);

  // Check surface indices
  EXPECT_EQ(gTTPtr->getSurfID(0), 1);
  EXPECT_EQ(gTTPtr->getSurfID(1), DAGMC::VOID_INDEX);

  // Check surfaces linked to volumes
  const std::vector<int> surfs = gTTPtr->getSurfs(0);
  ASSERT_EQ(surfs.size(), 1);
  EXPECT_EQ(surfs.at(0), 1);

  // Both vols should point at same surf
  const std::vector<int> surfs2 = gTTPtr->getSurfs(1);
  ASSERT_EQ(surfs2.size(), 1);
  EXPECT_EQ(surfs2.at(0), 1);

  // This volume doesn't exist: throw exception
  ASSERT_THROW(gTTPtr->getSurfs(2), std::out_of_range);

  // Check vols linked to surfaces
  const std::pair<int, int> volPairPtr = gTTPtr->getVolPair(0);
  EXPECT_EQ(volPairPtr.first, DAGMC::IMPLICIT_COMPLEMENT);
  EXPECT_EQ(volPairPtr.second, 1);

  // This surface doesn't exist: throw exception
  ASSERT_THROW(gTTPtr->getVolPair(1), std::out_of_range);

  //Done!

}

