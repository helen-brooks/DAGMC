#include "libmesh_test.hpp"
#include "GeomTopoTool.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class GeomTopoToolLMTest : public libMeshFileTest {
 protected:

  GeomTopoToolLMTest() :
    meshIsLoaded(false),
    filename("sphere.e")
  {};

  // Initalize variables for each test
  virtual void SetUp() override {
    libMeshFileTest::SetUp();
    meshIsLoaded = Read(filename);
    gTTPtr = std::make_shared<DAGMC::GeomTopoToolLM>(meshPtr);
  };

  // Data members required for tests
  bool meshIsLoaded;
  std::string filename;
  std::shared_ptr<DAGMC::GeomTopoTool> gTTPtr;

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
  EXPECT_TRUE(meshIsLoaded);

}

// Test that we can set up geometry
TEST_F(GeomTopoToolLMTest, setupGeometry) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);

  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_NE(gTTPtr, nullptr);
  ASSERT_TRUE(meshIsLoaded);

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
  const std::vector<intmax_t> surfs = gTTPtr->getSurfs(0);
  ASSERT_EQ(surfs.size(), 1);
  EXPECT_EQ(surfs.at(0), 1);

  // Both vols should point at same surf
  const std::vector<intmax_t> surfs2 = gTTPtr->getSurfs(1);
  ASSERT_EQ(surfs2.size(), 1);
  EXPECT_EQ(surfs2.at(0), 1);

  // This volume doesn't exist: throw exception
  ASSERT_THROW(gTTPtr->getSurfs(2), std::out_of_range);

  // Check vols linked to surfaces
  const DAGMC::SurfaceSenses& senses = gTTPtr->getVolPair(0);
  EXPECT_EQ(senses.forwards, DAGMC::IMPLICIT_COMPLEMENT);
  EXPECT_EQ(senses.backwards, 1);

  // This surface doesn't exist: throw exception
  ASSERT_THROW(gTTPtr->getVolPair(1), std::out_of_range);

  //Done!

}

