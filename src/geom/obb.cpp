#include "obb.hpp"

namespace DAGMC {

void  OrientedBoundingBox::construct_obb(const_element_iterator elemBegin,
                                         const_element_iterator elemEnd) {

  // Get all element vertices and put in matrix
  Matrix points = getPointsMatrix(elemBegin, elemEnd);

  // Construct covariance matrix and find eigenbasis
  Matrix basis;
  constructBasis(points, basis);

  // Find the extremal points along the basis vecs
  Vector minPoint;
  Vector maxPoint;
  findExtremalPoints(points, basis, minPoint, maxPoint);

  // Construct a box with the min/max point aligned with basis vecs
  constructBox(basis, minPoint, maxPoint);

}

// Construct a matrix of points
Matrix OrientedBoundingBox::getPointsMatrix(const_element_iterator elemBegin,
                                            const_element_iterator elemEnd) {

  // Declare 1X3 Matrix
  Matrix pointsMatrix(1, 3);

  // Loop over elements
  for (; elemBegin != elemEnd; ++elemBegin) {
    const Elem& element = **elemBegin;
    unsigned int nNodes = element.n_nodes();
    for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
      // Get the point associated with this node
      const Point& point = element.point(iNode);
      // Create a new row in matrix using initializer list.
      Matrix pointRow(1, 3) = {point(0), point(1), point(2)};
      // Insert row
      pointsMatrix.insert_rows(1, pointRow);
    }
  }
  return pointsMatrix;

}


// Construct a basis from a matrix of points
void OrientedBoundingBox::constructBasis(Matrix& points, Matrix& basis) {

  // Calculate the covariance matrix for the points
  Matrix covMat = cov(points);

  if (covMat.n_rows != DIM || covMat.n_cols != DIM)
    return;

  // Find the eigenvectors of the covariance matrix
  Vector eigvals;
  eigs_sym(eigvals, eigvec, covMat, DIM);
  if (eigvec.n_rows != DIM || eigven.n_cols != DIM)
    return;

  // Ensure columns are normalised to unity
  basis = normalise(eigvec);

}

void OrientedBoundingBox::findExtremalPoints(Matrix& points, Matrix& basis,
                                             Vector& minPoint, Vector& maxPoint) {

  if (basis.n_rows != DIM || basis.n_cols != DIM)
    return;
  if (points.n_cols != DIM)
    return;

  // Placeholder for min and max coeffs in basis
  std::array<double, DIM> minCoeffs;
  std::array<double, DIM> maxCoeffs;

  // points has a row for every vertex
  unsigned int nPoints = points.n_rows;
  for (unsigned int ipoint = 0; ipoint < nPoints; ++ipoint) {
    RowVec row = points(ipoint);
    // Need to tranpose
    Vector pvec = row.trans();
    // Loop over every basis vector
    for (unsigned int idim = 0; idim < DIM; ++idim) {
      Vector bvec = basis(idim);
      // Project out coeff along this basis vector
      double coeff = dot(pvec, bvec);

      if (ipoint == 0) {
        minCoeffs[idim] = coeff;
        maxCoeffs[idim] = coeff;
      } else {
        // Is this the min max coeff so far?
        if (coeff < minCoeffs[idim]) {
          minCoeff[idim] = coeff;
        }
        if (coeff > maxCoeffs[idim]) {
          maxCoeff[idim] = coeff;
        }
      }
    }// End loop over basis
  }// End loop over points


  // Now fetch the extremal points wrt the original axies
  minPoint.zeros();
  maxPoint.zeros();
  for (unsigned int idim = 0; idim < DIM; ++idim) {
    minPoint += minCoef[idim] * basis(idim);
    maxPoint += maxCoef[idim] * basis(idim);
  }


}


void OrientedBoundingBox::constructBox(Matrix& basis, Vector& minPoint, Vector& maxPoint) {

  // Find the vector to get between the extremal points
  Vector diff = maxPoint - minPoint;

  // Project the diff along the basis to get lengths of sides of box
  std::array<double, DIM> sides;
  for (unsigned int idim = 0; idim < DIM; ++idim) {
    sides[idim] = dot(diff, basis(idim));
  }

  std::array<Vector, 8> boxPoints;
  boxPoints[0] = minPoint;
  boxPoints[1] = boxPoints[0] + sides[0] * basis(0);
  boxPoints[2] = boxPoints[0] + sides[1] * basis(1);
  boxPoints[3] = boxPoints[1] + sides[1] * basis(1);
  boxPoints[4] = boxPoints[0] + sides[2] * basis(2);
  boxPoints[5] = boxPoints[4] + sides[0] * basis(0);
  boxPoints[6] = boxPoints[4] + sides[1] * basis(1);
  boxPoints[7] = boxPoints[5] + sides[1] * basis(1);

  //To-do check boxPoints[7] = maxPoint

}

}
