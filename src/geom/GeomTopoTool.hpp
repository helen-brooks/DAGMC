// C++ include files that we need
#include <iostream>
// // Functions to initialize the library.
// #include "libmesh/libmesh.h"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
// For node_iterator
#include "libmesh/mesh_base.h"
// For node
#include "libmesh/node.h"
// For elem
#include "libmesh/elem.h"

namespace DAGMC {

class GeomTopoTool {

 public:
  GeomTopoTool() {};
  ~GeomTopoTool() {};

  virtual bool setupGeometry() = 0;

  virtual void print() = 0;

 private:


};


class GeomTopoToolLM : public GeomTopoTool {

 public:


  GeomTopoToolLM(std::shared_ptr<libMesh::Mesh> meshPtrIn) :
    meshPtr(meshPtrIn), IMPLICIT_COMPLEMENT(-1) {
    std::cout << "Passed a pointer to a mesh." << std::endl;
  };

  ~GeomTopoToolLM() {};


  bool setupGeometry() override;

  void print() override;

 private:

  // Shared pointer to the mesh data
  std::shared_ptr<libMesh::Mesh> meshPtr;

  // Booking for topological relationship heirarchy
  std::map<int, std::vector<int> > vol_to_surfs;
  std::map<int, std::pair<int, int> > surf_to_vols;
  const int IMPLICIT_COMPLEMENT;

};

}// End namespace DAGMC
