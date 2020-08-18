#include <iostream>
#include <memory>
#include <gtest/gtest.h>
#include "obbtree.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class OBBTreeTest : public ::testing::Test {

 protected:

  OBBTreeTest() : libMeshException(false) {
    tol = 0.000000001;
  };

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
  // Tolerance for double comparisons
  double tol;
  // Save whether libMesh raised an exception
  bool libMeshException;

};

//---------------------------------------------------------------------------//

class SimpleOBBTreeTest : public OBBTreeTest {
 protected:

  SimpleOBBTreeTest() {
    nFaces = 2;
    nNodes = 4;
    nNodesPerFace = 3;
    points.resize(nNodes);
    conn.resize(nFaces);
  };

  // Initalize variables for each test
  virtual void SetUp() override {
    OBBTreeTest::SetUp();
    initMesh();
  };

  //Override
  virtual void initMesh() {
    // Need to generate a simple mesh.
    // We'll use a tetrahedron
    if (meshPtr != nullptr) {

      // Clear any prior data
      meshPtr->clear();

      libmesh_try {
        // Set dim
        meshPtr->set_mesh_dimension(2);
        meshPtr->set_spatial_dimension(3);

        // Reserve space
        meshPtr->reserve_elem(nFaces);
        meshPtr->reserve_nodes(nNodes);

        // Add points
        // Designed such that mean vec is (0,0,0)
        points.at(0) = libMesh::Point(-0.5, -0.5, 0.0);
        points.at(1) = libMesh::Point(-0.5, 0.5, 0.0);
        points.at(2) = libMesh::Point(0.5, 0.5, 0.0);
        points.at(3) = libMesh::Point(0.5, -0.5, 0.0);
        for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
          meshPtr->add_point(points.at(iNode), iNode);
        }

        // Add elems
        conn.at(0) = {0, 1, 2};
        conn.at(1) = {2, 3, 0};
        for (unsigned int iElem = 0; iElem < nFaces; iElem++) {
          libMesh::Elem* elemptr = meshPtr->add_elem(libMesh::Elem::build_with_id(libMesh::TRI3, iElem));
          for (unsigned int iNode = 0; iNode < nNodesPerFace; iNode++) {
            elemptr->set_node(iNode) = meshPtr->node_ptr(conn[iElem][iNode]);
          }
        }

        meshPtr->prepare_for_use();
      }
      libmesh_catch(libMesh::LogicError & e) {
        libMeshException = true;
      }
    }
  };

  // Data members
  unsigned int nFaces;
  unsigned int nNodes;
  unsigned int nNodesPerFace;
  std::vector<libMesh::Point> points;
  std::vector< std::vector< unsigned int> > conn;
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

