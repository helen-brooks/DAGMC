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

  // Find out how many degerate sides there are
  // 0 : box
  // 1 : square
  // 2 : line segment
  // 3 : point
  unsigned int nDegenerate() const {return degenAxes.size(); };

  // Get the ith degenerate direction. Corresponds to
  // col index in basis
  unsigned int getDegenDir(unsigned int i) const {
    return degenAxes.at(i);
  }

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
  // Face, edges and corners (or endpoints for degenerate boxes)
  // are inclusive
  // N.B. Will return false if box not sane
  bool containsPoint(Vector point) const;

  // Return the ith basis vector
  // Throw std::logic_error if out of bounds
  const Vector basisVec(unsigned int i) const {
    return basis.col(i);
  }

  // Define a status code for sanity checks
  enum BoxStatus { success, faildim, failnorm, failorth,
                   failordered
                 };

 private:

  // On construction, perform some sanity checks
  BoxStatus getStatus(Matrix& M, Vector& min, Vector& max,
                      std::vector<unsigned int>& degen) const;

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
  // Save if box is degenerate along any axes
  std::vector<unsigned int> degenAxes;

};


}
