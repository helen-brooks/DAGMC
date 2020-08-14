#include "box.hpp"
#include "container.hpp"

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
    method(methodIn) {
    elemsPtr = std::make_shared<DAGMC::ElemConstItContainer>(elemBegin, elemEnd);
    construct_obb();
  }

  bool isConstructed() { return (box != nullptr); }
  bool isSane() { return (box == nullptr) ? false : box->isSane(); };
  int  status() { return (box == nullptr) ? Box::failunknown : box->getBoxStatus(); };

  // Return a const reference to box
  const Box& getBox() { return *box; };

  // Geom queries
  bool containsPoint(const Vector& point) {
    return (box == nullptr) ? false : box->containsPoint(point);
  };

  // Overloaded version. Can also be used for libMesh::Node
  // since that inherits from libMesh::Point
  bool containsPoint(const libMesh::Point& pointLM) {
    Vector pointVec = { pointLM(0), pointLM(1), pointLM(2) };
    return containsPoint(pointVec);
  };


 private:
  // Construct an OBB for an element set
  void construct_obb();

  const ConstructMethod method;

  // The actual box
  std::shared_ptr<Box> box;

  // Storage for elements
  std::shared_ptr<ElemContainer> elemsPtr;

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
void constructBasisDiscrete(ElemContainer& elems,
                            Matrix& basis, Matrix& points);

// Construct basis using the covariance with a continuous average
// over each element (assumes tris).
// This is the second method presented in "OBBTree:
// A Hierarchical Structure for Rapid Interference Detection",
// Gottschalk, Lin, and Manocha
// section 4
void constructBasisCont(ElemContainer& elems,
                        Matrix& basis, Matrix& points);

// Construct a matrix of points
void getPointsMatrix(ElemContainer& elems, Matrix& points);

// Compute all statistics needed to calculate covariance from a
// set of elems. Assumes tris.
void getElemStats(ElemContainer& elems,
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
