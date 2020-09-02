#include "libmesh_interface.hpp"
#include "io_factory.hpp"

namespace DAGMC {

// *****************************************************************************
// PUBLIC METHODS
// *****************************************************************************

InternalMesh::InternalMesh(int argc, const char* const* argv) {
  // Initialize libMesh and handle any exceptions
  libmesh_try {
    // Initialize the library
    initPtr = std::make_shared<libMesh::LibMeshInit>(argc, argv);
    // Create the mesh
    if (initPtr != nullptr) {
      meshPtr = std::make_shared<libMesh::Mesh>(initPtr->comm());
    }
  }
  libmesh_catch(libMesh::LogicError & e) {
    return;
  }
};

LibMeshInterface::LibMeshInterface(int argc, const char* const* argv) {
  container = std::make_shared<InternalMesh>(argc, argv);
}

LibMeshInterface::LibMeshInterface(libMesh::MeshBase& meshRefIn) {
  container = std::make_shared<ExternalMesh>(meshRefIn);
}


bool LibMeshInterface::load(std::string filename) {

  if (!loadMesh(filename)) {
    std::cerr << "Failed to load mesh." << std::endl;
    return false;
  } else
    return loadSenseData(filename);
}

// *****************************************************************************
// PRIVATE METHODS
// *****************************************************************************

bool LibMeshInterface::loadMesh(std::string filename) {

  bool meshLoaded = false;
  libmesh_try {
    // Read in the mesh data
    if (!meshIsNull()) {
      // Clear any prior data
      mesh().clear();
      mesh().read(filename);
      meshLoaded = mesh().is_prepared();
    }
  }
  libmesh_catch(libMesh::LogicError & e) {
    return false;
  }
  return meshLoaded;
}

bool LibMeshInterface::loadSenseData(std::string filename) {

  // Create a reader
  std::shared_ptr<IOBase> reader = getIOPtr(filename);
  if (reader == nullptr)
    return false;

  // Open file for reading
  if (!reader->open(filename, IOBase::READ)) {
    std::cerr << "Could not open " << filename << " for reading." << std::endl;
    return false;
  }

  if (!reader->readAttributes(attributes)) {
    std::cerr << "Failed to read attributes from file " << filename << std::endl;
    return false;
  } else {
    return true;
  }

}


}
