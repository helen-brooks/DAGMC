#ifndef DAG_MESH_INTERFACE_HPP
#define DAG_MESH_INTERFACE_HPP

#include "Common.hpp"
#include "Types.hpp"

namespace DAGMC {

// Structures for Mesh Atrributes

// Struct for the volumes for which this surf has forwards/backwards sense
struct SurfaceSenses {
  intmax_t forwards;
  intmax_t backwards;
};

// Intended to be a generic extensible way of storing
// mesh attribute data that is not stored in the mesh itself.
struct MeshAttributes {
  std::map<intmax_t, SurfaceSenses > senseData;
  // Add more here as required
};

// Wrapper to agnosticize whether mesh is internal or external
template < class MeshType >
class MeshContainer {
 public:
  MeshContainer<MeshType>() {};
  ~MeshContainer<MeshType>() {};
  virtual MeshType& mesh() = 0;
  virtual const MeshType& const_mesh() const = 0;
  virtual MeshType* ptr() { return nullptr; };
  virtual std::shared_ptr<MeshType> sptr() {
    throw std::runtime_error("Mesh instance is not defined as a shared pointer !");
    return nullptr;
  }
  virtual bool isNull() { return false; };

};


class MeshInterfaceBase {

 public:
  MeshInterfaceBase() {};
  ~MeshInterfaceBase() {};

  // Read and write to file
  virtual bool load(std::string filename) = 0;
  virtual bool write(std::string filename) = 0;

  // Finish the setup of the geometry from an open file
  virtual bool finish_setup() = 0;

  // Find the geometry sets
  virtual bool setup_geom() = 0;
  virtual bool setup_indices() = 0;

  // Return the error code of the last operation
  virtual ErrorCode code() = 0;

};

template < class MeshType >
class MeshInterface : public MeshInterfaceBase {

 public:

  MeshInterface<MeshType>() {};
  ~MeshInterface<MeshType>() {};

  bool meshIsNull() {
    if (container == nullptr)
      return true;
    else
      return container->isNull();
  };

  MeshType& mesh() {
    return container->mesh();
  };
  const MeshType& const_mesh() {
    return container->const_mesh();
  };

  MeshType* mesh_ptr() { return container->ptr(); };
  std::shared_ptr<MeshType> mesh_sptr() { return container->sptr(); };

 protected:

  // Container for the mesh
  std::shared_ptr<MeshContainer<MeshType> > container;
};

};


#endif
