#include "obb.hpp"
#include <limits>

namespace DAGMC {

//---------------------------------------------------------------------------//
// OrientedBoundingBox class methods
//---------------------------------------------------------------------------//

bool OrientedBoundingBox::containsElem(libMesh::dof_id_type id) const {

  // Fetch a local iterator
  std::shared_ptr<ElemIterator> it = elemsPtr->getIterator();

  bool found = false;
  // Loop over elements
  const libMesh::Elem* elemptr;
  while (it->getNext(elemptr)) {
    if (elemptr->id() == id) {
      found = true;
      break;
    }
  }

  return found;
}


void  OrientedBoundingBox::construct_obb() {

  if (elemsPtr == nullptr)
    return;
  //Fetch a reference to container
  const ElemContainer& elems = *elemsPtr;

  // Find basis vectors for the box
  Matrix basis;
  Matrix points;
  constructBasis(elems, basis, points);

  // Find the extremal points along the basis vecs
  Vector minPoint;
  Vector maxPoint;
  OBBUtils::findExtremalPoints(points, basis, minPoint, maxPoint);

  // Construct a box with the min/max point aligned with basis vecs
  box = std::make_shared<Box>(maxPoint, minPoint, basis);

}

void OrientedBoundingBox::constructBasis(const ElemContainer& elems,
                                         Matrix& basis,
                                         Matrix& points) {

  if (method == ConstructMethod::CONT) {
    OBBUtils::constructBasisCont(elems, basis, points, meanPoint);
  } else if (method == ConstructMethod::DISCRETE) {
    OBBUtils::constructBasisDiscrete(elems, basis, points, meanPoint);
  } else {
    //Just get points and an axis-aligned basis
    OBBUtils::getPointsMatrix(elems, points, meanPoint);
    basis.eye(DIM, DIM);
  }

}

bool OrientedBoundingBox::setChildren() {

  // Check construction was successful
  if (!isSane())
    return false;

  // Partition elements
  std::vector<std::shared_ptr<ElemContainer> > partitions;
  getPartitions(partitions);

  // If partitions were found, now create
  // child OBBs for each partition
  for (auto& part : partitions) {
    // Shouldn't happen
    if (part == nullptr)
      return false;

    std::shared_ptr<TreeNode> child = getChild(part);

    // Save
    if (child != nullptr)
      children.push_back(child);
    // error
    else
      return false;

  }

  // Done
  return true;
}

std::shared_ptr<TreeNode> OrientedBoundingBox::getChild(std::shared_ptr<ElemContainer> elemsPtr) {

  std::shared_ptr<TreeNode> child;
  try {
    child =  std::make_shared<OrientedBoundingBox>(elemsPtr, method, shared_from_this());
  } catch (std::bad_weak_ptr& e) {
    // This can happen if OBB was not created through make_shared
    std::cerr << "Error: " << e.what() << " in " << __func__ << std::endl;
    child = nullptr;
  }
  return child;

}

void OrientedBoundingBox::getPartitions(std::vector<std::shared_ptr<ElemContainer> >& partitions) const {

  // Sort basis by box side length
  std::vector<unsigned int> order;
  box->getBasisOrder(order);

  // No non-degnerate dirs, or something went wrong
  if (order.empty())
    return;

  // Fetch a local iterator
  std::shared_ptr<ElemIterator> it = elemsPtr->getIterator();

  // Attempt to subdivide perpendicular to sides of decreasing length
  for (auto& ivec : order) {

    // Splitting plane is defined as having normal of basisVector,
    // and contains mean point vector.
    const Vector bVec = box->basisVec(ivec);

    //Calculate the projection of the mean point vector on this side
    double meanCoord = arma::dot(meanPoint, bVec);

    // Partition elements according to which side of the
    // splitting plane they are
    std::set<const libMesh::Elem*> elemsPlus;
    std::set<const libMesh::Elem*> elemsMinus;

    // Send iterator back to start
    it->reset();

    // Loop over elements
    const libMesh::Elem* elemptr;
    while (it->getNext(elemptr)) {
      //Get element midpoint
      Vector midpoint = getElemMidpoint(elemptr);

      //Find projection along basis vector relative to mean
      double coord = arma::dot(midpoint, bVec) - meanCoord;

      // Sort (inclusive: elems in plane go into elemsMinus)
      if (coord > 0.) {
        elemsPlus.insert(elemptr);
      } else {
        elemsMinus.insert(elemptr);
      }
    }

    // If either set is empty, box is indivisible along this dir.
    if (elemsPlus.empty() || elemsMinus.empty())
      continue;
    else {

      // Create containers for these elements
      std::shared_ptr<ElemContainer> contPlus
        = std::make_shared<ElemConstPtrContainer>(elemsPlus);
      std::shared_ptr<ElemContainer> contMinus
        = std::make_shared<ElemConstPtrContainer>(elemsMinus);

      //Save
      partitions.push_back(contPlus);
      partitions.push_back(contMinus);

      // Done. Exit loop
      break;
    }
  }
}

Vector OrientedBoundingBox::getElemMidpoint(const libMesh::Elem* elemPtr) const {

  Vector midpoint;
  Matrix points;
  if (method == ConstructMethod::CONT) {
    double area;
    OBBUtils::getSingleElemStats(elemPtr, points, midpoint, area);
  } else {
    OBBUtils::getSingleElemPoints(elemPtr, points, midpoint);
  }
  return midpoint;
}


//---------------------------------------------------------------------------//
// OBBUtils namespace methods
//---------------------------------------------------------------------------//

void OBBUtils::constructBasisCont(const ElemContainer& elems,
                                  Matrix& basis,
                                  Matrix& points,
                                  Vector& mean) {

  // Compute statistics on this element set
  std::vector<double> areas;
  OBBUtils::getElemStats(elems, areas, mean, points);

  // Construct the covariance matrix from the statistics
  Matrix cov;
  OBBUtils::calcCov(areas, mean, points, cov);

  // Compute basis
  OBBUtils::constructBasisFromCov(cov, basis);

}

void OBBUtils::constructBasisCont(const ElemContainer& elems,
                                  Matrix& basis,
                                  Matrix& points) {
  Vector dummy;
  OBBUtils::constructBasisCont(elems, basis, points, dummy);
}


// Construct basis using the covariance of the nodes of elements
void OBBUtils::constructBasisDiscrete(const ElemContainer& elems,
                                      Matrix& basis,
                                      Matrix& points,
                                      Vector& mean) {

  // Get all element vertices and put in matrix
  OBBUtils::getPointsMatrix(elems, points, mean);

  // Construct the covariance matrix for the points
  // Pass transpose: function expects entries are rows
  Matrix covMat = arma::cov(points.t());

  // Compute basis
  OBBUtils::constructBasisFromCov(covMat, basis);

}

void OBBUtils::constructBasisDiscrete(const ElemContainer& elems,
                                      Matrix& basis,
                                      Matrix& points) {
  Vector dummy;
  constructBasisDiscrete(elems, basis, points, dummy);
}

// Construct a matrix of points from the nodes of a set of elements
void OBBUtils::getPointsMatrix(const ElemContainer& elems,
                               Matrix& points,
                               Vector& mean) {

  // Fetch a local iterator
  std::shared_ptr<ElemIterator> it = elems.getIterator();
  points.reset();
  mean.zeros(DIM);

  libmesh_try {

    // Loop over elements
    const libMesh::Elem* elemPtr;
    unsigned int nElems = 0;
    while (it->getNext(elemPtr)) {
      nElems++;

      // Insert points
      Vector elemMean;
      OBBUtils::getSingleElemPoints(elemPtr, points, elemMean);

      // Add contribution to mean
      mean += elemMean;

    }
    // Correctly normalise mean vector to number of elems
    mean *= (1.0 / double(nElems));

  }
  libmesh_catch(libMesh::LogicError & e) {
    points.reset();
  }
}
void OBBUtils::getPointsMatrix(const ElemContainer& elems,
                               Matrix& points) {
  Vector dummy;
  getPointsMatrix(elems, points, dummy);
}

void OBBUtils::getSingleElemPoints(const libMesh::Elem* elemPtr,
                                   Matrix& points,
                                   Vector& elemMean) {
  elemMean.zeros(DIM);
  if (elemPtr == nullptr)
    return;

  unsigned int nNodes = elemPtr->n_nodes();
  for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
    // Get the point associated with this node
    const libMesh::Point& point = elemPtr->point(iNode);
    // Create a new vector using initializer list.
    Vector pvec = {point(0), point(1), point(2)};
    // Insert column
    if (points.n_cols == 0) {
      points = pvec;
    } else {
      points.insert_cols(points.n_cols, pvec);
    }
    // Add contribution to mean
    elemMean += pvec;
  }
  //Normalise mean
  elemMean *= 1.0 / double(nNodes);
}

void OBBUtils::getElemStats(const ElemContainer& elems,
                            std::vector<double>& areas,
                            Vector& mean,
                            Matrix& points) {
  // Fetch a local iterator
  std::shared_ptr<ElemIterator> it = elems.getIterator();
  points.reset();
  areas.clear();
  mean.zeros(DIM);

  libmesh_try {

    // Loop over elements
    const libMesh::Elem* elemPtr;
    while (it->getNext(elemPtr)) {

      // Fetch data for this element, and append points
      double area;
      Vector elemMean;
      OBBUtils::getSingleElemStats(elemPtr, points, elemMean, area);

      // Add contribution to mean
      mean += elemMean;

      // Save area
      areas.push_back(area);

    } // End loop over elems

    // Fetch number of elems
    unsigned int nElems = areas.size();

    // Correctly normalise mean vector
    mean *= (1.0 / double(nElems));

  }
  libmesh_catch(libMesh::LogicError & e) {
    points.reset();
    mean.reset();
    areas.clear();
  }

}

void OBBUtils::getSingleElemStats(const libMesh::Elem* elemPtr,
                                  Matrix& points,
                                  Vector& elemMean,
                                  double& area) {

  elemMean.zeros(DIM);

  if (elemPtr == nullptr)
    return;

  unsigned int nNodes = elemPtr->n_nodes();
  // Require elments are tris
  libMesh::ElemType type = elemPtr->type();
  if (type != libMesh::ElemType::TRI3 || nNodes != 3) {
    points.reset();
    return;
  }

  unsigned int offset = points.n_cols;

  // Get points
  for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
    // Get the point associated with this node
    const libMesh::Point& point = elemPtr->point(iNode);
    // Get point
    Vector pvec = {point(0), point(1), point(2)};
    // Insert column
    if (points.n_cols == 0) {
      points = pvec;
    } else {
      points.insert_cols(points.n_cols, pvec);
    }

    // Compute contribution to mean
    elemMean += pvec;
  }

  const Vector& p = points.col(offset + 0);
  const Vector& q = points.col(offset + 1);
  const Vector& r = points.col(offset + 2);

  // Compute area
  area = 0.5 * arma::norm(arma::cross(q - p, r - p));

  // Normalise mean
  elemMean /= (area * 6.0);
}

void OBBUtils::constructBasisFromCov(Matrix& cov, Matrix& basis) {
  // Check the covariance matrix is square and symmetric
  if (cov.n_rows != DIM || cov.n_cols != DIM)
    return;
  else if (!cov.is_symmetric())
    return;

  // Find the eigenvectors of the covariance matrix
  Vector eigvals;
  Matrix eigvec;
  bool solved = arma::eig_sym(eigvals, eigvec, cov);
  if (!solved)
    return;
  if (eigvec.n_rows != DIM || eigvec.n_cols != DIM)
    return;

  // Ensure columns are normalised to unity
  basis = normalise(eigvec);

  // Armadillo convention: eigvectors ordered by increasing eigenvalues
  // Prefer decreasing so degenerate dims are last
  basis = arma::fliplr(basis);

}

void OBBUtils::calcCov(std::vector<double>& areas, Vector& mean,
                       Matrix& points, Matrix& cov) {

  cov.reset();

  // Check dimensions are self-consistent
  unsigned int nElems = areas.size();
  if (mean.n_rows != DIM)
    return;
  if (points.n_rows != DIM || points.n_cols != 3 * nElems)
    return;

  // Initialise cov matrix
  cov.zeros(DIM, DIM);

  for (unsigned int iElem = 0; iElem < nElems; iElem++) {

    // Fetch elem points and subtract mean
    unsigned int iPoint = 3 * iElem;
    Vector pBar = points.col(iPoint) - mean;
    Vector qBar = points.col(iPoint + 1) - mean;
    Vector rBar = points.col(iPoint + 2) - mean;
    Vector sum = pBar + qBar + rBar;

    // Compute outer products
    Matrix outerSum = arma::kron(sum, sum.t());
    outerSum       += arma::kron(pBar, pBar.t());
    outerSum       += arma::kron(qBar, qBar.t());
    outerSum       += arma::kron(rBar, rBar.t());

    // Normalise
    outerSum *= 1.0 / (24.0 * double(nElems) * areas.at(iElem));

    // Add contribution to covariance
    cov +=  outerSum;

  }
  //Done
}


void OBBUtils::findExtremalPoints(Matrix& points, Matrix& basis,
                                  Vector& minPoint, Vector& maxPoint) {

  // Initialise min/max points
  minPoint.zeros(DIM);
  maxPoint.zeros(DIM);

  if (basis.n_rows != DIM || basis.n_cols != DIM)
    return;
  if (points.n_rows != DIM)
    return;

  // points has a col for every vertex
  unsigned int nPoints = points.n_cols;
  if (nPoints == 0)
    return;

  // Find the coeffs of the min/max points along the basis vectors
  double MAX_DOUBLE = std::numeric_limits<double>::max();
  std::vector<double> minCoeffs(DIM, MAX_DOUBLE);
  std::vector<double> maxCoeffs(DIM, -MAX_DOUBLE);

  for (unsigned int ipoint = 0; ipoint < nPoints; ++ipoint) {
    Vector pvec = points.col(ipoint);

    // Loop over every basis vector
    for (unsigned int idim = 0; idim < DIM; ++idim) {
      Vector bvec = basis.col(idim);
      // Project out coeff along this basis vector
      double coeff = dot(pvec, bvec);

      // Is this the min or max coeff so far?
      if (coeff < minCoeffs.at(idim)) {
        minCoeffs.at(idim) = coeff;
      }
      if (coeff > maxCoeffs.at(idim)) {
        maxCoeffs.at(idim) = coeff;
      }

    }// End loop over basis
  }// End loop over points

  // Now fetch the extremal points wrt the original axes
  for (unsigned int idim = 0; idim < DIM; ++idim) {
    minPoint += minCoeffs.at(idim) * basis.col(idim);
    maxPoint += maxCoeffs.at(idim) * basis.col(idim);
  }
}

}
