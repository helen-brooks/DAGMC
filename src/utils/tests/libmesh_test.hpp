#ifdef LIBMESH
#include <gtest/gtest.h>
#include "Libmesh.hpp"
#include <iostream>
#include <memory>

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class libMeshTest : public ::testing::Test {

 protected:

  libMeshTest() :
    libMeshException(false),
    tol(0.000000001) {};

  virtual void SetUp() override {
    initlibMesh();
  };
  virtual void TearDown() override {};

  void initlibMesh() {
    //Emulate dummy command line args since required by libmesh
    int argc_dummy = 1;
    char dummych[] = "dummy";
    char* argv_dummy[] = { dummych };

    // Reset pointer
    meshPtr.reset();

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

// Testing class where the mesh is read from file
class libMeshFileTest : public libMeshTest {
 protected:

  libMeshFileTest() {};

  bool Read(std::string filename) {
    bool meshLoaded = false;
    libmesh_try {
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
      return false;
    }
    return meshLoaded;
  }

};

// Testing class where the mesh is built manually
class libMeshSimpleTest : public libMeshTest {

 protected:

  // Initalize variables for each test
  virtual void SetUp() override {
    libMeshTest::SetUp();
    initMesh();
  };

  //Override
  virtual void initMesh() {

    // Generate a simple mesh.
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
        for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
          meshPtr->add_point(pointsLM.at(iNode), iNode);
        }

        // Add elems
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

  // Be sure to set pointsLM and conn
  std::vector<libMesh::Point> pointsLM;
  std::vector< std::vector< unsigned int> > conn;
  unsigned int nFaces;
  unsigned int nNodes;
  unsigned int nNodesPerFace;

};
#endif
