#ifndef DAG_OBB_HPP
#define DAG_OBB_HPP

#ifdef LIBMESH

#include "Common.hpp"
#include "Tree.hpp"
#include "Box.hpp"
#include "Container.hpp"

namespace DAGMC {

// 3D space is assumed throughout
const unsigned int DIM = 3;

enum ConstructMethod { CONT = 0, DISCRETE, ALIGNED };

// Class to construct a box containing a set of libMesh elements
// Default is oriented, but can build aligned box
class OrientedBoundingBox : public TreeNode {

 public:
  OrientedBoundingBox(const_element_iterator elemBegin,
                      const_element_iterator elemEnd,
                      ConstructMethod methodIn = ConstructMethod::CONT,
                      std::shared_ptr<TreeNode> parentIn = nullptr) :
    method(methodIn),
    TreeNode(parentIn) {
    elemsPtr = std::make_shared<ElemConstItContainer>(elemBegin, elemEnd);
    construct_obb();
  }

  OrientedBoundingBox(std::shared_ptr<ElemContainer> elemsPtrIn,
                      ConstructMethod methodIn = ConstructMethod::CONT,
                      std::shared_ptr<TreeNode> parentIn = nullptr) :
    method(methodIn),
    elemsPtr(elemsPtrIn),
    TreeNode(parentIn) {
    construct_obb();
  };

  // Return information on success of construction
  bool isConstructed() const override { return (box != nullptr); }
  bool isSane() const { return (box == nullptr) ? false : box->isSane(); };
  int  status() const { return (box == nullptr) ? Box::failunknown : box->getBoxStatus(); };

  // Return a const reference to box
  const Box& getBox() const { return *box; };

  // Geom queries
  bool containsPoint(const Vector& point) const {
    return (box == nullptr) ? false : box->containsPoint(point);
  };
  // Overloaded version. Can also be used for libMesh::Node
  // since that inherits from libMesh::Point
  bool containsPoint(const libMesh::Point& pointLM) const {
    Vector pointVec = { pointLM(0), pointLM(1), pointLM(2) };
    return containsPoint(pointVec);
  };

  // Not a geometric query - might be that an element could fit inside box
  // - but instead checks if this elem lives inside elemcontainer
  bool containsElem(libMesh::dof_id_type id) const;

  // Get a local iterator for the container
  std::shared_ptr<ElemIterator> getIterator() const {
    if (elemsPtr == nullptr)
      return nullptr;
    else
      return elemsPtr->getIterator();
  }

 private:
  // Construct an OBB for an element set
  void construct_obb();

  // Construct basis for the OBB and internally set mean vector
  void constructBasis(const ElemContainer& elems, Matrix& basis, Matrix& points);

  // Partition elements if possible and create new child OBBs
  bool setChildren() override;

  // Create a child box for a set of elems
  std::shared_ptr<TreeNode> getChild(std::shared_ptr<ElemContainer> elemsPtr);

  // Implmement subdivision algorithm from "OBBTree:
  // A Hierarchical Structure for Rapid Interference Detection",
  // Gottschalk, Lin, and Manocha, section 4
  void getPartitions(std::vector<std::shared_ptr<ElemContainer> >& partitions) const;

  Vector getElemMidpoint(const libMesh::Elem* elemPtr) const;

  // Method for constructing box
  const ConstructMethod method;

  // Mean vector of points
  Vector meanPoint;

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
void constructBasisDiscrete(const ElemContainer& elems,
                            Matrix& basis, Matrix& points, Vector& mean);
void constructBasisDiscrete(const ElemContainer& elems,
                            Matrix& basis, Matrix& points);

// Construct basis using the covariance with a continuous average
// over each element (assumes tris).
// This is the second method presented in "OBBTree:
// A Hierarchical Structure for Rapid Interference Detection",
// Gottschalk, Lin, and Manocha
// section 4
void constructBasisCont(const ElemContainer& elems,
                        Matrix& basis, Matrix& points, Vector& mean);
void constructBasisCont(const ElemContainer& elems,
                        Matrix& basis, Matrix& points);

// Construct a matrix of points (optionally fetch mean)
void getPointsMatrix(const ElemContainer& elems, Matrix& points, Vector& mean);
void getPointsMatrix(const ElemContainer& elems, Matrix& points);

//  Insert points for a single element and calc mean point
void getSingleElemPoints(const libMesh::Elem* elemPtr,
                         Matrix& points,
                         Vector& elemMean);

// Compute all statistics needed to calculate covariance from a
// set of elems. Assumes tris.
void getElemStats(const ElemContainer& elems,
                  std::vector<double>& areas, Vector& mean,
                  Matrix& points);

// Insert points for single element into points, and fetch stats
// Assumes tris.
void getSingleElemStats(const libMesh::Elem* elemPtr,
                        Matrix& points,
                        Vector& elemMean,
                        double& area);

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

#endif
#endif
