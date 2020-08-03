#include <armadillo>

namespace DAGMC {

typedef armadillo::mat Matrix;
typedef armadillo::vec Vector;
typedef armadillo::rowvec RowVec;
const unsigned int DIM = 3;

class OrientedBoundingBox {

 public:
  OrientedBoundingBox(const_element_iterator elemBegin,
                      const_element_iterator elemEnd) {
    construct_obb(elemBegin, elemEnd);
  }


 private:
  // Construct an OBB for an element set
  void construct_obb(const_element_iterator elemBegin,
                     const_element_iterator elemEnd);

  // Construct a matrix of points
  Matrix getpointsMatrix(const_element_iterator elemBegin,
                         const_element_iterator elemEnd);

  // Construct a basis from a matrix of points
  void constructBasis(Matrix& points, Matrix& basis);

  // Given a set of points and a basis, find the most extreme coefficients
  // along each basis vector
  void findExtremalPoints(Matrix& points, Matrix& basis,
                          Vector& minPoint, Vector& maxPoint);

  // Given a basis, and a min/max point, construct a box
  void constructBox(Matrix& basis, Vector& minPoint, Vector& maxPoint);

}

}
