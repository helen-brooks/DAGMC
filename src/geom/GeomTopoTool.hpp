#ifndef DAG_GTT_HPP
#define DAG_GTT_HPP

#ifdef LIBMESH
#include "LibmeshInterface.hpp"

namespace DAGMC {

// Use maximum integer width so that we are guaranteed to fit libMesh ID types
const intmax_t IMPLICIT_COMPLEMENT = -1;
const intmax_t VOID_INDEX = INTMAX_MIN;

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
  virtual intmax_t getVolID(unsigned int index) = 0;
  virtual intmax_t getSurfID(unsigned int index) = 0;

  // Return constant reference to surfaces and volumes
  // associated with a given index
  virtual const std::vector<intmax_t>& getSurfs(unsigned int index) = 0;
  virtual const SurfaceSenses& getVolPair(unsigned int index) = 0;

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
  intmax_t getVolID(unsigned int index) override;
  intmax_t getSurfID(unsigned int index) override;

  // Return constant reference to surfaces and volumes
  // associated with a given index
  const std::vector<intmax_t>& getSurfs(unsigned int index) override;
  const SurfaceSenses& getVolPair(unsigned int index) override;

 private:

  // Shared pointer to the mesh data
  std::shared_ptr<libMesh::Mesh> meshPtr;

  // Booking for topological relationship heirarchy
  std::map<intmax_t, std::vector<intmax_t> > vol_to_surfs;
  std::map<intmax_t, SurfaceSenses > surf_to_vols;
};

}// End namespace DAGMC

#endif

#endif
