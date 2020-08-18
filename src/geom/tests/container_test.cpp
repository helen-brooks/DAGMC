#include "libmesh_test.hpp"
#include "container.hpp"

// Super-lightweight test for the container class.

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class ContainerTest : public libMeshFileTest {

 protected:

  ContainerTest() : meshIsLoaded(false), filename("cube.e") {};

  virtual void SetUp() override {
    libMeshFileTest::SetUp();
    meshIsLoaded = Read(filename);
  };

  bool meshIsLoaded;
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
  ASSERT_TRUE(meshIsLoaded);

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
