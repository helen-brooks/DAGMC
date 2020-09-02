#ifndef DAG_LIBMESH_INTERFACE_HPP
#define DAG_LIBMESH_INTERFACE_HPP

#include "libmesh.hpp"
#include "mesh_interface.hpp"

namespace DAGMC {

// Wrapper to agnosticize whether mesh is internal or external
class MeshContainer {
 public:
  MeshContainer() {};
  ~MeshContainer() {};
  virtual libMesh::MeshBase& mesh() = 0;
  virtual const libMesh::MeshBase& const_mesh() const = 0;
  virtual bool isNull() { return false; };
};

// Saves a reference to an external mesh
class ExternalMesh : public MeshContainer {

 public:
  ExternalMesh(libMesh::MeshBase& meshRefIn) : _mesh(meshRefIn) {};

  libMesh::MeshBase& mesh() override { return _mesh; };
  const libMesh::MeshBase& const_mesh() const override { return _mesh; };

 private:

  libMesh::MeshBase& _mesh;

};

// This container owns the mesh and calls its contructor
class InternalMesh : public MeshContainer {

 public:

  InternalMesh(int argc, const char* const* argv);

  libMesh::MeshBase& mesh() override { return *meshPtr; };
  const libMesh::MeshBase& const_mesh() const override { return *meshPtr; };
  bool isNull() override { return (meshPtr == nullptr); };

 private:

  // Pointers to libMesh objects
  std::shared_ptr<libMesh::LibMeshInit> initPtr;
  std::shared_ptr<libMesh::MeshBase> meshPtr;

};


class LibMeshInterface : public MeshInterface {
 public:

  // This constructor creates an internal mesh
  LibMeshInterface(int argc, const char* const* argv);

  // This constructor createa an external mesh
  LibMeshInterface(libMesh::MeshBase& meshRefIn);

  ~LibMeshInterface() {};

  // Load data from filename into memory
  virtual bool load(std::string filename) override;

  // Get a copy of sense data
  std::map<intmax_t, SurfaceSenses > getSenseData() {
    return attributes.senseData;
  };

  bool meshIsNull() {
    if (container == nullptr)
      return true;
    else
      return container->isNull();
  };

  libMesh::MeshBase& mesh() {
    return container->mesh();
  };
  const libMesh::MeshBase& const_mesh() {
    return container->const_mesh();
  };

 private:

  // Load the primary mesh data (elements, nodes, etc)
  bool loadMesh(std::string filename);

  // Load surface to volume sense relationships
  bool loadSenseData(std::string filename);

  // Container for the mesh
  std::shared_ptr<MeshContainer> container;

  // Additional mesh attributes
  MeshAttributes attributes;

};

}

#endif
