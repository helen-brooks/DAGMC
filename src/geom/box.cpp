#include "box.hpp"
#include <limits>
#include <algorithm>

namespace DAGMC {


// Generic constructor: oriented box
Box::Box(Vector max,  Vector min,  Matrix M) :
  dim(max.n_rows), maxPoint(max), minPoint(min), basis(M) {
  // Sanity checks
  status = getStatus(basis, minPoint, maxPoint, degenAxes);
  sane = (status == BoxStatus::success);
}

// Specialist constructor: aligned box
Box::Box(Vector max,  Vector min) :
  dim(max.n_rows), maxPoint(max), minPoint(min) {
  // Set basis to identity matrix
  basis.eye(dim, dim);
  // Sanity checks
  status = getStatus(basis, minPoint, maxPoint, degenAxes);
  sane = (status == BoxStatus::success);
}

// Specialist constructor: aligned box with
// min point on origin
Box::Box(Vector max) :
  dim(max.n_rows), maxPoint(max) {
  minPoint.zeros(dim);
  // Set basis to identity matrix
  basis.eye(dim, dim);
  // Sanity checks
  status = getStatus(basis, minPoint, maxPoint, degenAxes);
  sane = (status == BoxStatus::success);
}

Box::BoxStatus Box::getStatus(Matrix& M, Vector& min, Vector& max,
                              std::vector<unsigned int>& degen) const {

  // Check that all dimensions are consistent
  if (min.n_rows != dim || max.n_rows != dim || M.n_rows != dim
      || M.n_cols != dim)
    return BoxStatus::faildim;


  // Check that max > min in new basis
  Vector diff = max - min;

  // Check that our basis vectors are orthonormal
  // Allow some finite tolerance for numerical effects.
  double tol = 1e-9;
  for (int icol = 0; icol < dim; ++icol) {
    Vector veci = M.col(icol);
    // Check normalised
    if (fabs(norm(veci) - 1.0) > tol) {
      return BoxStatus::failnorm;
    }
    // Check orthogonal
    for (int jcol = icol + 1; jcol < dim; ++jcol) {
      Vector vecj = M.col(jcol);
      if (fabs(dot(veci, vecj)) > tol) {
        return BoxStatus::failorth;
      }
    }

    // Project out diff along this basis vec to get side
    double side = dot(diff, veci);

    // Check that points have correct ordering (max>min)
    if (side < 0)
      return BoxStatus::failordered;

    // Check if box would have a zero-length side
    // i.e. box is degnerate
    if (side == 0) {
      degen.push_back(icol);
    }

  }
  // Our box passed all tests.
  return BoxStatus::success;
}

// Query if ray intersects box
// See "An Introduction to Ray Tracing" book for more details
bool Box::intersectsRay(const Vector& orig, const Vector& dir,
                        double& tEnter, double& tExit) const {

  // Do nothing for broken box
  if (!sane)
    return false;

  // Ray is defined by R = orig + dir*t
  // Box faces are defined by point.n = 0

  // Initialise intersection points at +/- inf
  double tfar = std::numeric_limits<double>::max();
  double tnear = -std::numeric_limits<double>::max();

  // Loop over basis vectors
  for (int icol = 0; icol < dim; ++icol) {
    // Vector is face normal n
    Vector faceNorm = basisVec(icol);

    // Get all dot products we need
    double dirProj = dot(faceNorm, dir);
    double minface = dot(minPoint, faceNorm);
    double maxface = dot(maxPoint, faceNorm);
    double raystart = dot(orig, faceNorm);

    // Check if ray is parallel to box faces
    if (dirProj == 0.) {
      // We'll only hit box if the ray origin is between
      // the two faces with this normal
      if (raystart > maxface ||
          raystart < minface)
        return false;
    } else {
      // Not parallel to box faces with this normal
      // -> must intersect.
      // Calculate instersection points
      // t = ((point - orig).n)/(dir.n)
      double t1;
      double t2;
      if (dirProj > 0.) {
        t1 = (minface - raystart) / dirProj;
        t2 = (maxface - raystart) / dirProj;
      } else {
        t2 = (minface - raystart) / dirProj;
        t1 = (maxface - raystart) / dirProj;
      }

      // Compare intersection points to previous plane pair
      // Require ordering to be bounded by planes
      if (t1 > tfar)
        return false;
      if (t2 < tnear)
        return false;

      // Save pair for next comparison
      tfar  = std::min(tfar, t2);
      tnear = std::max(tnear, t1);
    }

  } //End for loop

  // If we got to here, plane intersection points
  // must be on box surface
  tEnter = tnear;
  tExit  = tfar;
  return true;
};

bool Box::intersectsRay(const Vector& orig, const Vector& dir) const {
  double tIn(0.);
  double tOut(0.);
  return intersectsRay(orig, dir, tIn, tOut);
}


// // Query if this box intersects another box
// bool Box::intersectsBox(const Box& box) const {

//   // Do nothing for broken box
//   if(!sane) return false;

// }

// Query if box contains point
bool Box::containsPoint(const Vector& point) const {

  // Do nothing for broken box
  if (!sane)
    return false;

  // Loop over basis vectors
  for (int icol = 0; icol < dim; ++icol) {

    // Basis vector is normal to pair of faces
    Vector faceNorm = basisVec(icol);

    // Get all dot products we need
    double pos     = dot(point,   faceNorm);
    double minface = dot(minPoint, faceNorm);
    double maxface = dot(maxPoint, faceNorm);

    //Is the point between this pair of faces?
    if (pos > maxface || pos < minface)
      return false;
  }

  // Got to here: point is bounded by all faces
  return true;
}

void Box::getBasisOrder(std::vector<unsigned int>& order) const {

  order.clear();
  // Save side length -> basis vector to sort automatically
  // vector is protection against degnerate sides
  std::map< double, std::vector< unsigned int> > sides;
  getSides(sides);

  // Loop over sides in decreasing side length
  for (auto side = sides.rbegin(); side != sides.rend(); ++side) {
    for (auto& ivec : side->second) {
      order.push_back(ivec);
    }
  }

  // Something went wrong
  if (order.size() != (dim - nDegenerate())) {
    order.clear();
  }

}

void Box::getSides(std::map< double, std::vector< unsigned int> >& sides) const {

  sides.clear();
  Vector diff = maxPoint - minPoint;
  for (unsigned int idim = 0; idim < dim; idim++) {
    double length = arma::dot(diff, basis.col(idim));
    // Ignore degenerate directions.
    if (length == 0.)
      continue;
    if (sides.find(length) == sides.end()) {
      sides[length] = std::vector<unsigned int>();
    }
    sides[length].push_back(idim);
  }
}

}
