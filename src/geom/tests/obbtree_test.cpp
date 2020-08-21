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

class SimpleSquareOBBTreeTest : public libMeshSimpleTest {
 protected:

  SimpleSquareOBBTreeTest() {
    // Square of four tris
    nFaces = 4;
    nNodes = 5;
    nNodesPerFace = 3;

    pointsLM.resize(nNodes);
    pointsLM.at(0) = libMesh::Point(0., 0., 0.);
    pointsLM.at(1) = libMesh::Point(1., 1., 0.);
    pointsLM.at(2) = libMesh::Point(1., -1., 0.);
    pointsLM.at(3) = libMesh::Point(-1., -1., 0.);
    pointsLM.at(4) = libMesh::Point(-1., 1., 0.);

    conn.resize(nFaces);
    conn.at(0) = {0, 1, 2};
    conn.at(1) = {0, 2, 3};
    conn.at(2) = {0, 3, 4};
    conn.at(3) = {0, 4, 1};
  };

};


class FileOBBTreeTest : public libMeshFileTest {

 protected:

  FileOBBTreeTest(std::string filenameIn) :
    meshIsLoaded(false),
    filename(filenameIn) {};
  ~FileOBBTreeTest() {};

  virtual void SetUp() override {
    libMeshFileTest::SetUp();
    meshIsLoaded = Read(filename);
  };

  bool meshIsLoaded;
  std::string filename;

};

class CubeOBBTreeTest : public FileOBBTreeTest {

 protected:
  CubeOBBTreeTest() : nFaces(6), nElemInFace(26), FileOBBTreeTest("cube.e") {};
  ~CubeOBBTreeTest() {};

  unsigned int nFaces;
  unsigned int nElemInFace;
};

class SphereOBBTreeTest : public FileOBBTreeTest {

 protected:
  SphereOBBTreeTest() : FileOBBTreeTest("sphere.e") {};
  ~SphereOBBTreeTest() {};

};

//---------------------------------------------------------------------------//
// EXTRA FUNCTIONS DECLARATIONS
//---------------------------------------------------------------------------//

// Want to use these with several fixtures, so declare as standalone

// Helper function to organise tree nodes by depth
bool flattenTree(DAGMC::OBBTree& tree,
                 std::map< unsigned int, std::vector< std::shared_ptr<DAGMC::TreeNode> > >& treemap);

void doTreeTest(DAGMC::const_element_iterator elBeg,
                DAGMC::const_element_iterator elEnd,
                unsigned int nFaces);

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

TEST_F(SimpleSquareOBBTreeTest, Constructor) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());
  EXPECT_EQ(meshPtr->n_elem(), nFaces);
  EXPECT_EQ(meshPtr->n_nodes(), nNodes);

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  // Do the test
  doTreeTest(elBeg, elEnd, nFaces);

}

TEST_F(CubeOBBTreeTest, Constructor) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  // Get surfaces
  std::set< libMesh::subdomain_id_type > surf_ids;
  meshPtr->subdomain_ids(surf_ids);
  EXPECT_EQ(surf_ids.size(), nFaces);

  // 2) Loop over surfaces
  for (auto isurf : surf_ids) {

    DAGMC::const_element_iterator elBeg
      = meshPtr->active_subdomain_elements_begin(isurf);
    DAGMC::const_element_iterator elEnd
      = meshPtr->active_subdomain_elements_end(isurf);

    // Do the test
    doTreeTest(elBeg, elEnd, nElemInFace);
  }

}

TEST_F(SphereOBBTreeTest, Constructor) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);
  ASSERT_TRUE(meshPtr->is_prepared());

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  // Do the test
  doTreeTest(elBeg, elEnd, meshPtr->n_elem());

}

//---------------------------------------------------------------------------//
// EXTRA FUNCTIONS: DEFINITIONS
//---------------------------------------------------------------------------//

// Return false if any nullptrs encountered
bool flattenTree(DAGMC::OBBTree& tree,
                 std::map< unsigned int, std::vector< std::shared_ptr<DAGMC::TreeNode> > >& treemap) {

  treemap.clear();

  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  if (root == nullptr)
    return false;

  unsigned int depth = 0;
  treemap[depth] = std::vector<std::shared_ptr<DAGMC::TreeNode> >(1, root);

  std::vector<std::shared_ptr<DAGMC::TreeNode> > children
    = root->getChildren();

  while (!children.empty()) {
    // Save children from previous level
    depth++;
    treemap[depth] = children;

    // Get next level children
    children.clear();
    for (auto& nodeptr : treemap[depth]) {
      if (nodeptr == nullptr) {
        std::cerr << "Null pointer encountered at depth: "
                  << nodeptr << std::endl;
        return false;
      }
      if (!nodeptr->isLeaf()) {
        std::vector<std::shared_ptr<DAGMC::TreeNode> > nodeChildren
          = nodeptr->getChildren();
        children.insert(children.end(), nodeChildren.begin(), nodeChildren.end());
      }
    }
  }

  return true;
}

void doTreeTest(DAGMC::const_element_iterator elBeg,
                DAGMC::const_element_iterator elEnd,
                unsigned int nFaces) {

  // Build a tree
  DAGMC::OBBTree tree(elBeg, elEnd);

  // Fetch root
  std::shared_ptr<DAGMC::TreeNode> root = tree.getRoot();
  ASSERT_NE(root, nullptr);

  std::map< unsigned int, std::vector< std::shared_ptr<DAGMC::TreeNode> > > treemap;
  // Construct a map of depth -> tree nodes
  // Check it worked
  ASSERT_TRUE(flattenTree(tree, treemap));
  ASSERT_FALSE(treemap.empty());
  ASSERT_TRUE(treemap.find(0) != treemap.end());
  ASSERT_FALSE(treemap.at(0).empty());
  EXPECT_EQ((treemap.at(0)).size(), 1);
  EXPECT_EQ((treemap.at(0)).at(0), root);

  // Now inspect tree starting from deepest nodes
  auto treeit = treemap.rbegin();
  auto rend = treemap.rend();
  rend--;
  unsigned int nLeaves = 0;
  std::set<libMesh::dof_id_type> elemIDs;
  for (; treeit != rend; ++treeit) {
    // Depth should be greater than zero
    EXPECT_GT(treeit->first, 0);
    EXPECT_FALSE((treeit->second).empty());
    for (auto& nodeptr : treeit->second) {
      EXPECT_FALSE(nodeptr->isRoot());
      if (nodeptr->isLeaf()) {
        nLeaves++;
        // Attempt to dynamically cast;
        std::shared_ptr<DAGMC::OrientedBoundingBox> obbptr
          = std::dynamic_pointer_cast<DAGMC::OrientedBoundingBox>(nodeptr);
        ASSERT_NE(obbptr, nullptr);
        std::shared_ptr<DAGMC::ElemIterator> it = obbptr->getIterator();

        // Leaf should contain exactly one elem
        const libMesh::Elem* elemptr;
        unsigned int nElems = 0;
        while (it->getNext(elemptr)) {
          nElems++;
          // Only expect one iteration
          EXPECT_EQ(nElems, 1);
          ASSERT_NE(elemptr, nullptr);
          elemIDs.insert(elemptr->id());
        }
      }
    }
  }


  EXPECT_EQ(nLeaves, nFaces);
  EXPECT_EQ(elemIDs.size(), nFaces);

  // Do a conventional loop over mesh elements to see ids match up
  DAGMC::const_element_iterator elIt = elBeg;
  for (; elIt != elEnd; elIt++) {
    libMesh::dof_id_type id = (**elIt).id();
    auto it = elemIDs.find(id);
    EXPECT_TRUE(it != elemIDs.end());
    // Don't expect to hit same id twice: delete
    if (it != elemIDs.end()) {
      elemIDs.erase(it);
    }
  }
  EXPECT_TRUE(elemIDs.empty());
}
