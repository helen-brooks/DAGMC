#include <gtest/gtest.h>
#include <iostream>
#include <memory>

#ifdef LIBMESH

#include "libmesh_interface.hpp"
#include "io_factory.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

class LibMeshInterfaceTest : public ::testing::Test {

 protected:

  LibMeshInterfaceTest() {
    filenames.push_back("sphere-withattributes.e");
  };

  void loadInternal() {
    //Emulate dummy command line args since required by libmesh
    int argc_dummy = 1;
    char dummych[] = "dummy";
    char* argv_dummy[] = { dummych };
    interface = std::make_shared<DAGMC::LibMeshInterface>(argc_dummy, argv_dummy);
  }

  void loadExternal(libMesh::MeshBase& meshRef) {
    interface = std::make_shared<DAGMC::LibMeshInterface>(meshRef);
  }

  // Data members required for tests
  std::shared_ptr<DAGMC::MeshInterface> interface;
  std::vector<std::string> filenames;

};

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS
//---------------------------------------------------------------------------//

TEST_F(LibMeshInterfaceTest, IOPtr) {

  std::shared_ptr<DAGMC::IOBase> ioptr =
      DAGMC::getIOPtr(filenames.at(0));

  EXPECT_NE(ioptr, nullptr);

}

TEST_F(LibMeshInterfaceTest, loadInternal) {

  loadInternal();
  ASSERT_NE(interface, nullptr);
  EXPECT_TRUE(interface->load(filenames.at(0)));

}

TEST_F(LibMeshInterfaceTest, loadExternal) {

  //Emulate dummy command line args since required by libmesh
  int argc_dummy = 1;
  char dummych[] = "dummy";
  char* argv_dummy[] = { dummych };

  // Initialize the library
  libMesh::LibMeshInit init(argc_dummy, argv_dummy);
  // Create the mesh
  libMesh::Mesh mesh(init.comm());

  loadExternal(mesh);
  ASSERT_NE(interface, nullptr);
  EXPECT_TRUE(interface->load(filenames.at(0)));

}

#endif
