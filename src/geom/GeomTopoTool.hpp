// C++ include files that we need
#include <iostream>
// Libmesh headers
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"

namespace DAGMC {

const int IMPLICIT_COMPLEMENT = -1;
const int VOID_INDEX = -999;

class GeomTopoTool {

 public:
  GeomTopoTool() {};
  ~GeomTopoTool() {};

  virtual bool setupGeometry() = 0;

  virtual void print() = 0;

  // Return number of entities
  virtual unsigned int nVols() = 0;
  virtual unsigned int nSurfs() = 0;

  // Convert from index to ID
  virtual int getVolID(unsigned int index) = 0;
  virtual int getSurfID(unsigned int index) = 0;

  // Return constant reference to surfaces and volumes
  // associated with a given index
  virtual const std::vector<int>& getSurfs(unsigned int index) = 0;
  virtual const std::pair<int, int>& getVolPair(unsigned int index) = 0;

 private:


};


class GeomTopoToolLM : public GeomTopoTool {

 public:


  GeomTopoToolLM(std::shared_ptr<libMesh::Mesh> meshPtrIn) :
    meshPtr(meshPtrIn) {};

  ~GeomTopoToolLM() {};

  bool setupGeometry() override;

  void print() override;

  // Note: answer includes implicit complement
  unsigned int nVols() override {
    return vol_to_surfs.size();
  };
  unsigned int nSurfs() override {
    return surf_to_vols.size();
  };

  // Convert from index to ID
  int getVolID(unsigned int index) override;
  int getSurfID(unsigned int index) override;

  // Return constant reference to surfaces and volumes
  // associated with a given index
  const std::vector<int>& getSurfs(unsigned int index) override;
  const std::pair<int, int>& getVolPair(unsigned int index) override;

 private:

  // Shared pointer to the mesh data
  std::shared_ptr<libMesh::Mesh> meshPtr;

  // Booking for topological relationship heirarchy
  std::map<int, std::vector<int> > vol_to_surfs;
  std::map<int, std::pair<int, int> > surf_to_vols;
  std::shared_ptr< std::vector<int> > dummy;

};

}// End namespace DAGMC
