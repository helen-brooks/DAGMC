#ifndef DAG_MESH_INTERFACE_HPP
#define DAG_MESH_INTERFACE_HPP

#include <memory>
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <cstdint>

namespace DAGMC {
class MeshInterface {

 public:

  MeshInterface() {};
  ~MeshInterface() {};

  virtual bool load(std::string filename) = 0;
  virtual bool write(std::string filename) = 0;

};

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
}

// Wrapper to agnosticize whether mesh is internal or external
template < class MeshType >
class MeshContainer {
 public:
  MeshContainer<MeshType>() {};
  ~MeshContainer<MeshType>() {};
  virtual MeshType& mesh() = 0;
  virtual const MeshType& const_mesh() const = 0;
  virtual bool isNull() { return false; };
};

#endif
