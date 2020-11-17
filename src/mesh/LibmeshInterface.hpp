#ifndef DAG_LIBMESH_INTERFACE_HPP
#define DAG_LIBMESH_INTERFACE_HPP

#include "MeshInterface.hpp"

#ifdef LIBMESH
#include "Libmesh.hpp"

namespace DAGMC {

// Saves a reference to an external mesh
class ExternalMesh : public MeshContainer<libMesh::MeshBase> {

 public:
  ExternalMesh(libMesh::MeshBase& meshRefIn) : _mesh(meshRefIn) {};

  libMesh::MeshBase& mesh() override { return _mesh; };
  const libMesh::MeshBase& const_mesh() const override { return _mesh; };
  libMesh::MeshBase* ptr() override { return &_mesh; };

 private:

  libMesh::MeshBase& _mesh;

};

// This container owns the mesh and calls its contructor
class InternalMesh : public MeshContainer<libMesh::MeshBase> {

 public:

  InternalMesh(int argc, const char* const* argv);

  libMesh::MeshBase& mesh() override { return *meshPtr; };
  const libMesh::MeshBase& const_mesh() const override { return *meshPtr; };
  libMesh::MeshBase* ptr() override { return meshPtr.get(); };
  std::shared_ptr<libMesh::MeshBase> sptr() override { return meshPtr; };
  bool isNull() override { return (meshPtr == nullptr); };

 private:

  // Pointers to libMesh objects
  std::shared_ptr<libMesh::LibMeshInit> initPtr;
  std::shared_ptr<libMesh::MeshBase> meshPtr;

};

class LibMeshInterface : public MeshInterface<libMesh::MeshBase> {
 public:

  // This constructor creates an internal mesh
  LibMeshInterface(int argc, const char* const* argv);

  // This constructor createa an external mesh
  LibMeshInterface(libMesh::MeshBase& meshRefIn);

  ~LibMeshInterface() {};

  // Load data from filename into memory
  bool load(std::string filename) override;
  bool write(std::string filename) override;

  // Finish the setup of the geometry from an open file
  bool finish_setup() override;

  // Find the geometry sets
  bool setup_geom() override;
  bool setup_indices() override;

  // Return the error code of the last operation
  ErrorCode code() override { return DAG_NOT_IMPLEMENTED; }

  // Get a copy of sense data
  std::map<intmax_t, SurfaceSenses > getSenseData() {
    return attributes.senseData;
  };

 private:

  // Load the primary mesh data (elements, nodes, etc)
  bool loadMesh(std::string filename);

  // Load surface to volume sense relationships
  bool loadSenseData(std::string filename);

  // Additional mesh attributes
  MeshAttributes attributes;

};

}
#endif

#endif