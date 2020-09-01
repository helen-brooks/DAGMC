#ifndef DAG_LIBMESH_INTERFACE_HPP
#define DAG_LIBMESH_INTERFACE_HPP

#include "libmesh.hpp"
#include "mesh_interface.hpp"

namespace DAGMC {

class LibMeshInterface : public MeshInterface {
 public:

  LibMeshInterface(int argc, const char* const* argv);
  ~LibMeshInterface() {};

  // Load data from filename into memory
  virtual bool load(std::string filename) override;

  // Get a copy of sense data
  std::map<intmax_t, SurfaceSenses > getSenseData() {
    return attributes.senseData;
  };

 private:

  // Load the primary mesh data (elements, nodes, etc)
  bool loadMesh(std::string filename);

  // Load surface to volume sense relationships
  bool loadSenseData(std::string filename);

  // Pointers to libMesh objects
  std::shared_ptr<libMesh::LibMeshInit> initPtr;
  std::shared_ptr<libMesh::Mesh> meshPtr;

  // Additional mesh attributes
  MeshAttributes attributes;

};

}

#endif
