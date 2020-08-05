#include <armadillo>

namespace DAGMC {

typedef arma::mat Matrix;
typedef arma::vec Vector;
typedef arma::rowvec RowVec;

// Class to perform geometrical queries on a generic box
class Box {

 public:
  // Generic constructor: oriented box
  Box(Vector max,  Vector min,  Matrix M);
  // Specialist constructor: aligned box
  Box(Vector max,  Vector min);
  // Specialist constructor: aligned box with
  // min point on origin
  Box(Vector max);

  // Default destructor
  ~Box() {};

  // Check construction succeeded
  bool isSane() const { return  sane; }
  int getBoxStatus() const { return int(status); };

  // Implements the ray - box intersection algorithm
  // from "An Introduction to Ray Tracing", contribution
  // "Essential Ray Tracing Algorithms" by Eric Haines,
  // section 4.
  // N.B. Will return false if box not sane
  // Optionally returns the distances to enter and exit the box
  // If box is behind ray, will return true, and tEnter/tExit
  // will be negative
  bool intersectsRay(const Vector& orig, const Vector& dir) const;
  bool intersectsRay(const Vector& orig, const Vector& dir,
                     double& tEnter, double& tExit) const;

  // // Query if this box intersects another box
  // bool intersectsBox(const Box& box) const;

  // Query if box contains point
  // N.B. Will return false if box not sane
  bool containsPoint(Vector point) const;

  // Return the ith basis vector
  // Throw std::logic_error if out of bounds
  const Vector basisVec(unsigned int i) const {
    return basis.col(i);
  }

  // Define a status code for sanity checks
  enum BoxStatus { success, faildim, failnorm, failorth,
                   faildegenerate, failordered
                 };

 private:

  // On construction, perform some sanity checks
  BoxStatus getStatus(Matrix& M, Vector& min, Vector& max) const;

  // Dimension of the box
  const unsigned int dim;

  // Extremal points of box with coefficients in the
  // original problem basis
  Vector maxPoint;
  Vector minPoint;

  // Matrix whose columns are the basis vectors
  // with coefficients in the original problem basis
  Matrix basis;

  // Save result of sanity checks
  BoxStatus status;
  bool sane;

};


}
