#include <gtest/gtest.h>
#include "container.hpp"

// Super-lightweight test for the container class.

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class ContainerTest : public ::testing::Test {

 protected:

  ContainerTest() : libMeshException(false),
    meshLoaded(false),
    filename("cube.e") {};

  virtual void SetUp() override {
    initlibMesh();
  };
  virtual void TearDown() override {};

  void initlibMesh() {
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
      // Read in the mesh data
      if (meshPtr != nullptr) {
        // Clear any prior data
        meshPtr->clear();
        meshPtr->read(filename);
        meshLoaded = meshPtr->is_prepared();
      }
    }
    libmesh_catch(libMesh::LogicError & e) {
      libMeshException = true;
      return;
    }
  }

  // Data members required for tests
  // Pointers to libMesh objects
  std::shared_ptr<libMesh::LibMeshInit> initPtr;
  std::shared_ptr<libMesh::Mesh> meshPtr;
  // Save whether libMesh raised an exception
  bool libMeshException;
  bool meshLoaded;
  std::string filename;
  std::shared_ptr<DAGMC::ElemContainer> dagElems;

};


//---------------------------------------------------------------------------//
// FIXTURE BASED TESTS
//---------------------------------------------------------------------------//

TEST_F(ContainerTest, Constructor) {

  ASSERT_FALSE(libMeshException);
  ASSERT_NE(initPtr, nullptr);
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshLoaded);

  // Get element iterators
  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();
  DAGMC::const_element_iterator elIt = elBeg;

  std::set<const libMesh::Elem*> elemptrs;
  std::vector<const libMesh::Elem*> elemptrsOrd;
  unsigned int nElems = meshPtr->n_elem();

  // Normal loop to set pointers for comparison
  for (; elIt != elEnd; ++elIt) {
    const libMesh::Elem& elem = **elIt;
    elemptrsOrd.push_back(&elem);
    elemptrs.insert(&elem);
  }
  // Need to check that setup didn't fail
  ASSERT_EQ(elemptrs.size(), nElems);
  ASSERT_EQ(elemptrsOrd.size(), nElems);

  // First constructor call
  dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elBeg, elEnd);

  // Get an iterator
  std::shared_ptr<DAGMC::ElemIterator> it = dagElems->getIterator();

  EXPECT_EQ(it->first(), &** elBeg);
  EXPECT_EQ(it->last(), &** elEnd);

  //Two attempts to check reset call works
  unsigned int attempts = 0;
  std::stringstream err;
  do {
    attempts++;
    err.str("");
    err << "Failed on attempt " << attempts;

    const libMesh::Elem* elemPtr;
    unsigned int count = 0;
    while (it->getNext(elemPtr)) {
      ASSERT_LT(count, nElems) << err.str();
      // Expect to hit the same elements in same order
      EXPECT_EQ(elemPtr, elemptrsOrd.at(count)) << err.str();
      count++;
    }
    EXPECT_EQ(count, nElems) << err.str();

    it->reset();

  } while (attempts < 2);

  // Second constructor call
  dagElems =
      std::make_shared<DAGMC::ElemConstPtrContainer>(elemptrs)
      ;
  // Get an iterator
  it = dagElems->getIterator();

  EXPECT_EQ(it->first(), *elemptrs.begin());
  EXPECT_EQ(it->last(), *elemptrs.end());

  //Two attempts to check reset call works
  attempts = 0;
  do {
    attempts++;
    err.str("");
    err << "Failed on attempt " << attempts;

    const libMesh::Elem* elemPtr;
    auto itCheck = elemptrs.begin();
    while (it->getNext(elemPtr)) {
      ASSERT_NE(itCheck, elemptrs.end());
      EXPECT_EQ(elemPtr, *itCheck);
      itCheck++;
    }

    //check reset call works
    it->reset();

  } while (attempts < 2);

}
