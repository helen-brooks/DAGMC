#include "obb.hpp"
#include <limits>

namespace DAGMC {

void  OrientedBoundingBox::construct_obb() {

  if (elemsPtr == nullptr)
    return;
  //Fetch a reference to container
  ElemContainer& elems = *elemsPtr;

  // Find basis vectors for the box
  Matrix basis;
  Matrix points;
  if (method == ConstructMethod::CONT) {
    OBBUtils::constructBasisCont(elems, basis, points);
  } else if (method == ConstructMethod::DISCRETE) {
    OBBUtils::constructBasisDiscrete(elems, basis, points);
  } else {
    //Just get points and an axis-aligned basis
    OBBUtils::getPointsMatrix(elems, points);
    basis.eye(DIM, DIM);
  }

  // Find the extremal points along the basis vecs
  Vector minPoint;
  Vector maxPoint;
  OBBUtils::findExtremalPoints(points, basis, minPoint, maxPoint);

  // Construct a box with the min/max point aligned with basis vecs
  box = std::make_shared<Box>(maxPoint, minPoint, basis);

}

void OBBUtils::constructBasisCont(ElemContainer& elems,
                                  Matrix& basis, Matrix& points) {

  // Compute statistics on this element set
  std::vector<double> areas;
  Vector mean;
  OBBUtils::getElemStats(elems, areas, mean, points);

  // Construct the covariance matrix from the statistics
  Matrix cov;
  OBBUtils::calcCov(areas, mean, points, cov);

  // Compute basis
  OBBUtils::constructBasisFromCov(cov, basis);

}


// Construct basis using the covariance of the nodes of elements
void OBBUtils::constructBasisDiscrete(ElemContainer& elems,
                                      Matrix& basis, Matrix& points) {

  // Get all element vertices and put in matrix
  OBBUtils::getPointsMatrix(elems, points);

  // Construct the covariance matrix for the points
  // Pass transpose: function expects entries are rows
  Matrix covMat = arma::cov(points.t());

  // Compute basis
  OBBUtils::constructBasisFromCov(covMat, basis);

}
// Construct a matrix of points from the nodes of a set of elements
void OBBUtils::getPointsMatrix(ElemContainer& elems,
                               Matrix& points) {
  elems.reset();
  points.reset();

  libmesh_try {
    // Loop over elements
    const libMesh::Elem* elemPtr;
    while (elems.getNext(elemPtr)) {
      unsigned int nNodes = elemPtr->n_nodes();
      for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
        // Get the point associated with this node
        const libMesh::Point& point = elemPtr->point(iNode);
        // Create a new row in matrix using initializer list.
        Vector pvec = {point(0), point(1), point(2)};
        // Insert column
        if (points.n_cols == 0) {
          points = pvec;
        } else {
          points.insert_cols(points.n_cols, pvec);
        }
      }
    }
  }
  libmesh_catch(libMesh::LogicError & e) {
    points.reset();
  }
}

void OBBUtils::getElemStats(ElemContainer& elems,
                            std::vector<double>& areas, Vector& mean,
                            Matrix& points) {
  elems.reset();
  points.reset();
  areas.clear();
  mean.zeros(DIM);

  libmesh_try {

    // Loop over elements
    const libMesh::Elem* elemPtr;
    while (elems.getNext(elemPtr)) {
      unsigned int nNodes = elemPtr->n_nodes();

      // Require elments are tris
      libMesh::ElemType type = elemPtr->type();
      if (type != libMesh::ElemType::TRI3 || nNodes != 3) {
        points.reset();
        return;
      }

      // Fetch coords of this element
      std::vector< Vector > elemPoints;
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
        elemPoints.push_back(pvec);
      }

      const Vector& p = elemPoints.at(0);
      const Vector& q = elemPoints.at(1);
      const Vector& r = elemPoints.at(2);

      // Compute area
      double area = 0.5 * arma::norm(arma::cross(q - p, r - p));
      areas.push_back(area);

      // Compute contribution to mean
      mean += (p + q + r) / (area * 6.0);

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
