#include "box.hpp"
// Libmesh headers
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"

namespace DAGMC {

const unsigned int DIM = 3;

typedef libMesh::MeshBase::const_element_iterator const_element_iterator;

// Class to construct a box containing a set of libMesh elements
class OrientedBoundingBox {

 public:
  OrientedBoundingBox(const_element_iterator elemBegin,
                      const_element_iterator elemEnd):
    elBegin(elemBegin), elEnd(elemBegin) {
    construct_obb(elemBegin, elemEnd);
  }

  bool isConstructed() { return (box != nullptr); }
  bool isSane() { return (box == nullptr) ? false : box->isSane(); };
  int  status() { return (box == nullptr) ? Box::failunknown : box->getBoxStatus(); };

  // Return a const reference to box
  const Box& getBox() { return *box; };

 private:
  // Construct an OBB for an element set
  void construct_obb(const_element_iterator elemBegin,
                     const_element_iterator elemEnd);

  // The actual box
  std::shared_ptr<Box> box;

  // Iterators to elements contained inside the box
  const_element_iterator elBegin;
  const_element_iterator elEnd;

  // A list of element IDs contained in the box?

};

// Some  helper functions
// Keep this separate so they can be used externally in principle
namespace OBBUtils {
// Construct a matrix of points
void getPointsMatrix(const_element_iterator elemBegin,
                     const_element_iterator elemEnd,
                     Matrix& points);

// Construct a basis from a matrix of points
void constructBasis(Matrix& points, Matrix& basis);

// Given a set of points and a basis, find the most extreme coefficients
// along each basis vector
void findExtremalPoints(Matrix& points, Matrix& basis,
                        Vector& minPoint, Vector& maxPoint);

};

}
