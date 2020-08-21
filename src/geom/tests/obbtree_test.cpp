#include "libmesh_test.hpp"
#include "obbtree.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class SimpleOBBTreeTest : public libMeshSimpleTest {
 protected:

  SimpleOBBTreeTest() {
    // Square of two tris
    nFaces = 2;
    nNodes = 4;
    nNodesPerFace = 3;

    pointsLM.resize(nNodes);
    pointsLM.at(0) = libMesh::Point(-0.5, -0.5, 0.0);
    pointsLM.at(1) = libMesh::Point(-0.5, 0.5, 0.0);
    pointsLM.at(2) = libMesh::Point(0.5, 0.5, 0.0);
    pointsLM.at(3) = libMesh::Point(0.5, -0.5, 0.0);

    conn.resize(nFaces);
    conn.at(0) = {0, 1, 2};
    conn.at(1) = {2, 3, 0};
  };

};

//---------------------------------------------------------------------------//
// FIXTURE BASED TESTS
//---------------------------------------------------------------------------//

TEST_F(SimpleOBBTreeTest, Constructor) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  // Build a tree and fetch root
  DAGMC::OBBTree tree(elBeg, elEnd);
  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  ASSERT_NE(root, nullptr);

  // number of uses:
  // original, 2xchildren, this copy
  EXPECT_EQ(root.use_count(), 4);

  EXPECT_TRUE(root->isRoot());
  EXPECT_FALSE(root->isLeaf());

  std::vector<std::shared_ptr<DAGMC::TreeNode> > children
    = root->getChildren();
  EXPECT_EQ(children.size(), 2);

  std::vector<unsigned int> elemBoxCount;
  elemBoxCount.resize(nFaces);
  for (auto& childptr : children) {
    // number of uses:
    // original, this copy
    EXPECT_EQ(childptr.use_count(), 2);

    EXPECT_FALSE(childptr->isRoot());
    EXPECT_TRUE(childptr->isLeaf());

    // Attempt to dynamically cast;
    std::shared_ptr<DAGMC::OrientedBoundingBox> obbptr
      = std::dynamic_pointer_cast<DAGMC::OrientedBoundingBox>(childptr);
    ASSERT_NE(obbptr, nullptr);

    // Both children should contain zero
    DAGMC::Vector zeroVec = {0., 0., 0.};
    EXPECT_TRUE(obbptr->containsPoint(zeroVec));

    // Neither child should contain 1,1,1
    DAGMC::Vector onesVec = {1., 1., 1.};
    EXPECT_FALSE(obbptr->containsPoint(onesVec));

    // Loop over elems
    bool boxContainsBoth = true;
    unsigned int nElemsInBox = 0;
    for (unsigned int iElem = 0; iElem < nFaces; iElem++) {
      bool boxContainsThis = obbptr->containsElem(iElem);
      if (!boxContainsThis) {
        boxContainsBoth = false;
      } else {
        elemBoxCount[iElem]++;
        nElemsInBox++;
      }
    }

    // Box should contain one and only one element
    EXPECT_FALSE(boxContainsBoth);
    EXPECT_EQ(nElemsInBox, 1);

  }

  //Each element should appear in one box only
  for (unsigned int iElem = 0; iElem < nFaces; iElem++) {
    EXPECT_EQ(elemBoxCount[iElem], 1);
  }

}


// Test an empty tree: repeated iterator
TEST_F(SimpleOBBTreeTest, DummyRepeated) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();

  // Build a tree and fetch root
  DAGMC::OBBTree tree(elBeg, elBeg);
  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  ASSERT_NE(root, nullptr);

  // Technically valid, just empty
  EXPECT_FALSE(root->isConstructed());

}

// Test an empty tree: inverted iterator
TEST_F(SimpleOBBTreeTest, DummyReversed) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  // Build a tree and fetch root
  DAGMC::OBBTree tree(elEnd, elBeg);
  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  ASSERT_NE(root, nullptr);

  EXPECT_FALSE(root->isConstructed());

}

// Test an empty tree: unordered iterators
TEST_F(SimpleOBBTreeTest, DummyUnordered) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elIt = elBeg;
  elIt++;

  // Build a tree and fetch root
  DAGMC::OBBTree tree(elIt, elBeg);
  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  ASSERT_NE(root, nullptr);

  EXPECT_FALSE(root->isConstructed());

}

// Test an empty tree: unmatched type of iterators
TEST_F(SimpleOBBTreeTest, DummyUnmatched) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->active_elements_end();

  // Build a tree and fetch root
  DAGMC::OBBTree tree(elBeg, elEnd);
  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  ASSERT_NE(root, nullptr);

  EXPECT_FALSE(root->isConstructed());
}


// Test a single node tree
TEST_F(SimpleOBBTreeTest, Leaf) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elIt = elBeg;
  elIt++;

  // Build a tree and fetch root
  DAGMC::OBBTree tree(elBeg, elIt);
  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  ASSERT_NE(root, nullptr);

  EXPECT_TRUE(root->isConstructed());
  EXPECT_TRUE(root->isRoot());
  EXPECT_TRUE(root->isLeaf());
}
