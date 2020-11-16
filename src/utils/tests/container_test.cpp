#ifdef LIBMESH

#include "libmesh_test.hpp"
#include "Container.hpp"

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

};

//---------------------------------------------------------------------------//
// FIXTURE BASED TESTS
//---------------------------------------------------------------------------//

TEST_F(ContainerTest, goodConstructor) {

  ASSERT_FALSE(libMeshException);
  ASSERT_NE(initPtr, nullptr);
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshIsLoaded);

  // Get element iterators
  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  std::shared_ptr<DAGMC::ElemContainer> dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elBeg, elEnd);
  EXPECT_TRUE(dagElems->isValid());

}

TEST_F(ContainerTest, getNext) {

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
  for (; elIt != elEnd; elIt++) {
    const libMesh::Elem& elem = **elIt;
    elemptrsOrd.push_back(&elem);
    elemptrs.insert(&elem);
  }
  // Need to check that setup didn't fail
  ASSERT_EQ(elemptrs.size(), nElems);
  ASSERT_EQ(elemptrsOrd.size(), nElems);

  //First constructor call
  std::shared_ptr<DAGMC::ElemContainer> dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elBeg, elEnd);
  EXPECT_TRUE(dagElems->isValid());

  // Get an iterator
  std::shared_ptr<DAGMC::ElemIterator> it
    = dagElems->getIterator();
  ASSERT_NE(it, nullptr);

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
      ASSERT_NE(elemPtr, nullptr);
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
      std::make_shared<DAGMC::ElemConstPtrContainer>(elemptrs);
  EXPECT_TRUE(dagElems->isValid());

  // Get an iterator
  it = dagElems->getIterator();
  ASSERT_NE(it, nullptr);

  //Two attempts to check reset call works
  attempts = 0;
  do {
    attempts++;
    err.str("");
    err << "Failed on attempt " << attempts;

    const libMesh::Elem* elemPtr;
    auto itCheck = elemptrs.begin();
    while (it->getNext(elemPtr)) {
      ASSERT_NE(itCheck, elemptrs.end()) << err.str();
      EXPECT_EQ(elemPtr, *itCheck) << err.str();
      itCheck++;
    }

    // check reset call works
    it->reset();

  } while (attempts < 2);

}

TEST_F(ContainerTest, BegIsEnd) {

  ASSERT_FALSE(libMeshException);
  ASSERT_NE(initPtr, nullptr);
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshIsLoaded);

  // Get element iterators
  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  // First constructor call
  std::shared_ptr<DAGMC::ElemContainer> dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elEnd, elBeg);
  EXPECT_FALSE(dagElems->isValid());

}

TEST_F(ContainerTest, EmptyContainer) {

  ASSERT_FALSE(libMeshException);
  ASSERT_NE(initPtr, nullptr);
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshIsLoaded);

  // Get element iterators
  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();

  // Repeated: should be valid but do nothing
  std::shared_ptr<DAGMC::ElemContainer> dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elBeg, elBeg);
  EXPECT_TRUE(dagElems->isValid());

  // Get an iterator
  std::shared_ptr<DAGMC::ElemIterator> it
    = dagElems->getIterator();
  ASSERT_NE(it, nullptr);

  const libMesh::Elem* elemPtr(nullptr);
  ASSERT_FALSE(it->getNext(elemPtr));
  EXPECT_EQ(elemPtr, nullptr);
}

TEST_F(ContainerTest, Unordered) {

  ASSERT_FALSE(libMeshException);
  ASSERT_NE(initPtr, nullptr);
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshIsLoaded);

  // Get element iterators
  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elIt = elBeg;
  elIt++;

  std::shared_ptr<DAGMC::ElemContainer> dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elIt, elBeg);
  ASSERT_FALSE(dagElems->isValid());

  // Get an iterator
  std::shared_ptr<DAGMC::ElemIterator> it
    = dagElems->getIterator();
  ASSERT_NE(it, nullptr);

  const libMesh::Elem* elemPtr;
  EXPECT_FALSE(it->getNext(elemPtr));
  EXPECT_EQ(elemPtr, nullptr);
}

TEST_F(ContainerTest, PredicateTypeTest) {

  ASSERT_FALSE(libMeshException);
  ASSERT_NE(initPtr, nullptr);
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshIsLoaded);

  // Get element iterators
  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  std::set< libMesh::subdomain_id_type > surf_ids;
  meshPtr->subdomain_ids(surf_ids);
  ASSERT_FALSE(surf_ids.empty());
  auto surfit = surf_ids.begin();
  libMesh::subdomain_id_type surf1 = *surfit;

  DAGMC::const_element_iterator elSurfBeg
    = meshPtr->active_subdomain_elements_begin(surf1);
  DAGMC::const_element_iterator elSurfEnd
    = meshPtr->active_subdomain_elements_end(surf1);

  // Different types of predicate: not valid
  std::shared_ptr<DAGMC::ElemContainer> dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elBeg, elSurfEnd);
  ASSERT_FALSE(dagElems->isValid());

  std::shared_ptr<DAGMC::ElemIterator> it
    = dagElems->getIterator();
  ASSERT_NE(it, nullptr);
  const libMesh::Elem* elemPtr;
  EXPECT_FALSE(it->getNext(elemPtr));
  EXPECT_EQ(elemPtr, nullptr);

  dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elSurfBeg, elEnd);
  ASSERT_FALSE(dagElems->isValid());

  it = dagElems->getIterator();
  ASSERT_NE(it, nullptr);
  EXPECT_FALSE(it->getNext(elemPtr));
  EXPECT_EQ(elemPtr, nullptr);

}

TEST_F(ContainerTest, PredicateArgTest) {

  ASSERT_FALSE(libMeshException);
  ASSERT_NE(initPtr, nullptr);
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshIsLoaded);

  std::set< libMesh::subdomain_id_type > surf_ids;
  meshPtr->subdomain_ids(surf_ids);
  ASSERT_FALSE(surf_ids.empty());
  ASSERT_GT(surf_ids.size(), 1);

  auto surfit = surf_ids.begin();
  libMesh::subdomain_id_type surf1 = *surfit;
  surfit++;
  libMesh::subdomain_id_type surf2 = *surfit;

  // Get element iterators
  DAGMC::const_element_iterator elSurf1Beg
    = meshPtr->active_subdomain_elements_begin(surf1);
  DAGMC::const_element_iterator elSurf1End
    = meshPtr->active_subdomain_elements_end(surf1);
  DAGMC::const_element_iterator elSurf2Beg
    = meshPtr->active_subdomain_elements_begin(surf2);
  DAGMC::const_element_iterator elSurf2End
    = meshPtr->active_subdomain_elements_end(surf2);
  DAGMC::const_element_iterator elSurf2It
    = elSurf2Beg;
  elSurf2It++;

  // Same type of predicate, but different subdomain: not valid
  std::shared_ptr<DAGMC::ElemContainer> dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elSurf1Beg, elSurf2End);
  ASSERT_FALSE(dagElems->isValid());

  std::shared_ptr<DAGMC::ElemIterator> it
    = dagElems->getIterator();
  ASSERT_NE(it, nullptr);

  const libMesh::Elem* elemPtr;
  EXPECT_FALSE(it->getNext(elemPtr));
  EXPECT_EQ(elemPtr, nullptr);

  dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elSurf1Beg, elSurf2It);
  ASSERT_FALSE(dagElems->isValid());
  it = dagElems->getIterator();
  ASSERT_NE(it, nullptr);
  EXPECT_FALSE(it->getNext(elemPtr));
  EXPECT_EQ(elemPtr, nullptr);

  dagElems =
      std::make_shared<DAGMC::ElemConstItContainer>(elSurf2Beg, elSurf1End);
  ASSERT_FALSE(dagElems->isValid());
  it = dagElems->getIterator();
  ASSERT_NE(it, nullptr);
  EXPECT_FALSE(it->getNext(elemPtr));
  EXPECT_EQ(elemPtr, nullptr);

}

#endif


