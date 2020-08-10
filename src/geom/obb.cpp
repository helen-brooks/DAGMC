#include "obb.hpp"
#include <limits>

namespace DAGMC {

void  OrientedBoundingBox::construct_obb(const_element_iterator elemBegin,
                                         const_element_iterator elemEnd) {

  // Get all element vertices and put in matrix
  Matrix points;
  OBBUtils::getPointsMatrix(elemBegin, elemEnd, points);

  // Construct covariance matrix and find eigenbasis
  Matrix basis;
  OBBUtils::constructBasis(points, basis);

  // Find the extremal points along the basis vecs
  Vector minPoint;
  Vector maxPoint;
  OBBUtils::findExtremalPoints(points, basis, minPoint, maxPoint);

  // Construct a box with the min/max point aligned with basis vecs
  box = std::make_shared<Box>(maxPoint, minPoint, basis);

}

// Construct a matrix of points
void OBBUtils::getPointsMatrix(const_element_iterator elemBegin,
                               const_element_iterator elemEnd,
                               Matrix& points) {
  points.reset();

  libmesh_try {
    // Loop over elements
    for (; elemBegin != elemEnd; ++elemBegin) {
      const libMesh::Elem& element = **elemBegin;
      unsigned int nNodes = element.n_nodes();
      for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
        // Get the point associated with this node
        const libMesh::Point& point = element.point(iNode);
        // Create a new row in matrix using initializer list.
        Matrix pointRow = {point(0), point(1), point(2)};
        // Insert row
        if (points.n_rows == 0) {
          points = pointRow;
        } else {
          points.insert_rows(points.n_rows, pointRow);
        }
      }
    }
  }
  libmesh_catch(libMesh::LogicError & e) {
    points.reset();
  }
}


// Construct a basis from a matrix of points
void OBBUtils::constructBasis(Matrix& points, Matrix& basis) {

  // Calculate the covariance matrix for the points
  Matrix covMat = cov(points);

  // Check the covariance matrix is square and symmetric
  if (covMat.n_rows != DIM || covMat.n_cols != DIM)
    return;
  else if (!covMat.is_symmetric())
    return;

  // Find the eigenvectors of the covariance matrix
  Vector eigvals;
  Matrix eigvec;
  bool solved = arma::eig_sym(eigvals, eigvec, covMat);
  if (!solved)
    return;
  if (eigvec.n_rows != DIM || eigvec.n_cols != DIM)
    return;

  // Ensure columns are normalised to unity
  basis = normalise(eigvec);

}

void OBBUtils::findExtremalPoints(Matrix& points, Matrix& basis,
                                  Vector& minPoint, Vector& maxPoint) {

  // Initialise min/max points
  minPoint.zeros(DIM);
  maxPoint.zeros(DIM);

  if (basis.n_rows != DIM || basis.n_cols != DIM)
    return;
  if (points.n_cols != DIM)
    return;

  // points has a row for every vertex
  unsigned int nPoints = points.n_rows;
  if (nPoints == 0)
    return;

  // Find the coeffs of the min/max points along the basis vectors
  double MAX_DOUBLE = std::numeric_limits<double>::max();
  std::vector<double> minCoeffs(DIM, MAX_DOUBLE);
  std::vector<double> maxCoeffs(DIM, -MAX_DOUBLE);
  for (unsigned int ipoint = 0; ipoint < nPoints; ++ipoint) {
    RowVec row = points.row(ipoint);
    // Need to tranpose
    Vector pvec = row.t();
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
    minPoint += minCoeffs.at(idim) * basis(idim);
    maxPoint += maxCoeffs.at(idim) * basis(idim);
  }
}

}
