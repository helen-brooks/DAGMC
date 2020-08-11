#include "box.hpp"
// Libmesh headers
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"

namespace DAGMC {

// 3D space is assumed throughout
const unsigned int DIM = 3;

enum ConstructMethod { CONT = 0, DISCRETE, ALIGNED };

typedef libMesh::MeshBase::const_element_iterator const_element_iterator;

// Class to construct a box containing a set of libMesh elements
// Default is oriented, but can build aligned box
class OrientedBoundingBox {

 public:
  OrientedBoundingBox(const_element_iterator elemBegin,
                      const_element_iterator elemEnd,
                      ConstructMethod methodIn = ConstructMethod::CONT) :
    method(methodIn),
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

  const ConstructMethod method;

  // The actual box
  std::shared_ptr<Box> box;

  // Iterators to elements contained inside the box
  const_element_iterator elBegin;
  const_element_iterator elEnd;

};

//  Functionally a typedef, just restricts the constructor
class AxisAlignedBoundingBox: public OrientedBoundingBox {

 public:
  AxisAlignedBoundingBox(const_element_iterator elemBegin,
                         const_element_iterator elemEnd) :
    OrientedBoundingBox(elemBegin, elemEnd, ConstructMethod::ALIGNED) {};

  ~AxisAlignedBoundingBox() {};

};

// Some  helper functions
// Keep this separate so they can be used externally in principle
namespace OBBUtils {

// Construct basis using the covariance of the nodes of elements.
// This is the first method presented in "OBBTree:
// A Hierarchical Structure for Rapid Interference Detection",
// Gottschalk, Lin, and Manocha
// section 4
void constructBasisDiscrete(const_element_iterator elemBegin,
                            const_element_iterator elemEnd,
                            Matrix& basis, Matrix& points);

// Construct basis using the covariance with a continuous average
// over each element (assumes tris).
// This is the second method presented in "OBBTree:
// A Hierarchical Structure for Rapid Interference Detection",
// Gottschalk, Lin, and Manocha
// section 4
void constructBasisCont(const_element_iterator elemBegin,
                        const_element_iterator elemEnd,
                        Matrix& basis, Matrix& points);

// Construct a matrix of points
void getPointsMatrix(const_element_iterator elemBegin,
                     const_element_iterator elemEnd,
                     Matrix& points);

// Compute all statistics needed to calculate covariance from a
// set of elems. Assumes tris.
void getElemStats(const_element_iterator elemBegin,
                  const_element_iterator elemEnd,
                  std::vector<double>& areas, Vector& mean,
                  Matrix& points);

// Compute the covariance given various statistics
void calcCov(std::vector<double>& areas, Vector& mean,
             Matrix& points, Matrix& cov);

// Construct a basis from a covariance matrix
void constructBasisFromCov(Matrix& cov, Matrix& basis);

// Given a set of points and a basis, find the most extreme coefficients
// along each basis vector
void findExtremalPoints(Matrix& points, Matrix& basis,
                        Vector& minPoint, Vector& maxPoint);

};

}
