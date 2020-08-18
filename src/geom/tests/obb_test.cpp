#include "libmesh_test.hpp"
#include "obb.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//

// An external helper function
bool checkBasis(const DAGMC::Matrix& basis, double tol, std::stringstream& errmsg);

//---------------------------------------------------------------------------//

class OBBTetTest : public libMeshSimpleTest {
 protected:

  OBBTetTest() {
    nFaces = 4;
    nNodes = 4;
    nNodesPerFace = 3;
    pointsLM.resize(nNodes);
    conn.resize(nFaces);

    // Designed such that mean vec is (0,0,0)
    pointsLM.at(0) = libMesh::Point(-0.5, -sqrt(3.) / 6., -sqrt(6.) / 12.);
    pointsLM.at(1) = libMesh::Point(0.0, sqrt(3.) / 3., -sqrt(6.) / 12.);
    pointsLM.at(2) = libMesh::Point(0.5, -sqrt(3.) / 6., -sqrt(6.) / 12.);
    pointsLM.at(3) = libMesh::Point(0.0,         0.0,  sqrt(6.) / 4.);

    conn.at(0) = {0, 1, 2};
    conn.at(1) = {1, 3, 0};
    conn.at(2) = {2, 3, 1};
    conn.at(3) = {0, 3, 2};

  };

};

//---------------------------------------------------------------------------//

class OBBFileTest : public libMeshFileTest {
 protected:

  OBBFileTest() {
    files.resize(2);
    files.at(0) = "cube.e";
    files.at(1) = "cube-rotate.e";

    vals.push_back(4.5);
    vals.push_back(7.2955549577);
  };

  // Save some information about mesh to help check later if it loaded correctly
  bool getSurfsAndElems() {

    surfs2Elems.clear();
    surfIDs.clear();

    libmesh_try {
      // Fetch surface ids
      std::set< libMesh::subdomain_id_type > surf_ids;
      meshPtr->subdomain_ids(surf_ids);

      // Loop over surfaces
      for (auto isurf : surf_ids) {

        // Get element ids
        std::set< libMesh::dof_id_type > elems;
        DAGMC::const_element_iterator elIt
        = meshPtr->active_subdomain_elements_begin(isurf);
        DAGMC::const_element_iterator elEnd
        = meshPtr->active_subdomain_elements_end(isurf);
        unsigned int nElems = 0;
        for (; elIt != elEnd; ++elIt) {
          const libMesh::Elem& elem = **elIt;
          elems.insert(elem.id());
          nElems++;
        }
        // Check all ids unique
        if (nElems != elems.size())
          return false;

        // Save elems
        surfs2Elems[isurf] = elems;
        surfIDs.push_back(isurf);
      }

      if (surfs2Elems.size() != surf_ids.size())
        return false;
      else
        return true;
    }
    libmesh_catch(libMesh::LogicError & e) {
      libMeshException = true;
      return false;
    }
  }

  // These methods are only used to check mesh was correctly loaded.
  unsigned int nSurfs() {
    return surfIDs.size();
  };
  unsigned int nElemsInSurf(unsigned int i) {

    libMesh::subdomain_id_type surfID = surfIDs.at(i);
    auto surfIt = surfs2Elems.find(surfID);
    if (surfIt != surfs2Elems.end()) {
      return (surfIt->second).size();
    } else
      return 0;
  }
  unsigned int nodesPerElemInSurf(unsigned int i) {
    libMesh::subdomain_id_type surfID = surfIDs.at(i);
    auto surfIt = surfs2Elems.find(surfID);
    if (surfIt != surfs2Elems.end()) {
      const libMesh::dof_id_type elemID = *((surfIt->second).begin());
      const libMesh::Elem* elemPtr = meshPtr->elem_ptr(elemID);
      if (elemPtr != nullptr)
        return elemPtr->n_nodes();
      else
        return 0;
    } else
      return 0;
  }


  std::map< libMesh::subdomain_id_type, std::set< libMesh::dof_id_type >  > surfs2Elems;
  std::vector<libMesh::subdomain_id_type> surfIDs;
  std::vector<std::string> files;
  std::vector<double> vals;
};

//---------------------------------------------------------------------------//
// FIXTURE BASED TESTS
//---------------------------------------------------------------------------//

TEST_F(OBBTetTest, SetUp) {
  ASSERT_FALSE(libMeshException);
  ASSERT_NE(meshPtr, nullptr);
  EXPECT_TRUE(meshPtr->is_prepared());
}

TEST_F(OBBTetTest, OBBUtilsTet) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();
  DAGMC::ElemConstItContainer elems(elBeg, elEnd);

  DAGMC::Matrix points;
  ASSERT_NO_THROW(DAGMC::OBBUtils::getPointsMatrix(elems, points));
  ASSERT_EQ(points.n_cols, nFaces * nNodesPerFace);
  ASSERT_EQ(points.n_rows, 3);

  // Check the points match up
  for (unsigned int iElem = 0; iElem < nFaces; iElem++) {
    for (unsigned int iNode = 0; iNode < nNodesPerFace; iNode++) {
      unsigned int iPoint = (conn.at(iElem)).at(iNode);
      libMesh::Point p;
      ASSERT_NO_THROW(p = pointsLM.at(iPoint));
      const DAGMC::Vector& pvec = points.col(iElem * nNodesPerFace + iNode);
      for (unsigned int iRow = 0; iRow < 3; iRow++) {
        EXPECT_EQ(p(iRow), pvec(iRow));
      }
    }
  }

  // Find the min/max points along x,y,z
  DAGMC::Matrix basis = { {1., 0., 0}, {0., 1., 0.}, {0., 0., 1.} };
  DAGMC::Vector min;
  DAGMC::Vector max;
  ASSERT_NO_THROW(DAGMC::OBBUtils::findExtremalPoints(points, basis, min, max));
  ASSERT_EQ(min.n_rows, 3);
  ASSERT_EQ(max.n_rows, 3);
  DAGMC::Vector minTest = {(pointsLM.at(0))(0),
                           (pointsLM.at(0))(1),
                           (pointsLM.at(0))(2)
                          };
  DAGMC::Vector maxTest = {(pointsLM.at(2))(0),
                           (pointsLM.at(1))(1),
                           (pointsLM.at(3))(2)
                          };

  EXPECT_LT(fabs(min(0) - minTest(0)), tol);
  EXPECT_LT(fabs(min(1) - minTest(1)), tol);
  EXPECT_LT(fabs(min(2) - minTest(2)), tol);
  EXPECT_LT(fabs(max(0) - maxTest(0)), tol);
  EXPECT_LT(fabs(max(1) - maxTest(1)), tol);
  EXPECT_LT(fabs(max(2) - maxTest(2)), tol);

  DAGMC::Matrix covMat = arma::cov(points.t());
  EXPECT_EQ(covMat.n_rows, 3);
  EXPECT_EQ(covMat.n_cols, 3);
  EXPECT_TRUE(covMat.is_symmetric());

  basis.reset();
  ASSERT_NO_THROW(DAGMC::OBBUtils::constructBasisFromCov(covMat, basis));
  std::stringstream errmsg;
  bool orthonormal = checkBasis(basis, tol, errmsg);
  EXPECT_TRUE(orthonormal) << errmsg.str();

  points.reset();
  basis.reset();
  covMat.reset();

  // Find element statistics
  std::vector<double> areas;
  DAGMC::Vector mean;
  ASSERT_NO_THROW(DAGMC::OBBUtils::getElemStats(elems, areas, mean, points));
  ASSERT_EQ(points.n_cols, nFaces * nNodesPerFace);
  ASSERT_EQ(points.n_rows, 3);
  ASSERT_EQ(mean.n_rows, 3);
  ASSERT_EQ(areas.size(), nFaces);

  // Check the points, areas and mean match up
  for (unsigned int iElem = 0; iElem < nFaces; iElem++) {
    for (unsigned int iNode = 0; iNode < nNodesPerFace; iNode++) {
      unsigned int iPoint = (conn.at(iElem)).at(iNode);
      libMesh::Point p;
      ASSERT_NO_THROW(p = pointsLM.at(iPoint));
      const DAGMC::Vector& pvec = points.col(iElem * nNodesPerFace + iNode);
      for (unsigned int iRow = 0; iRow < 3; iRow++) {
        EXPECT_EQ(p(iRow), pvec(iRow));
        EXPECT_LT(fabs(mean(iRow)), tol);
      }
    }
    double areaDiff = fabs(areas.at(iElem) - sqrt(3.) / 4.);
    EXPECT_LT(areaDiff, tol);
  }

  // Get covariance matrix from element stats
  covMat.reset();
  ASSERT_NO_THROW(DAGMC::OBBUtils::calcCov(areas, mean, points, covMat));

  // Get basis
  ASSERT_NO_THROW(DAGMC::OBBUtils::constructBasisFromCov(covMat, basis));
  errmsg.str("");
  orthonormal = checkBasis(basis, tol, errmsg);
  EXPECT_TRUE(orthonormal) << errmsg.str();

  // Compute basis directly through discrete method
  points.reset();
  basis.reset();
  errmsg.str("");
  ASSERT_NO_THROW(DAGMC::OBBUtils::constructBasisDiscrete(elems, basis, points));
  orthonormal = checkBasis(basis, tol, errmsg);
  EXPECT_TRUE(orthonormal) << errmsg.str();

  // Compute basis directly through discrete method
  points.reset();
  basis.reset();
  errmsg.str("");
  ASSERT_NO_THROW(DAGMC::OBBUtils::constructBasisCont(elems, basis, points));
  orthonormal = checkBasis(basis, tol, errmsg);
  EXPECT_TRUE(orthonormal) << errmsg.str();

  // Compute basis for a single face
  DAGMC::const_element_iterator el = elBeg;
  el++;
  elems = DAGMC::ElemConstItContainer(elBeg, el);
  basis.reset();
  errmsg.str("");
  ASSERT_NO_THROW(DAGMC::OBBUtils::constructBasisDiscrete(elems, basis, points));
  ASSERT_EQ(points.n_cols, nNodesPerFace);
  ASSERT_EQ(points.n_rows, 3);
  orthonormal = checkBasis(basis, tol, errmsg);
  EXPECT_TRUE(orthonormal) << errmsg.str();

  points.reset();
  basis.reset();
  errmsg.str("");
  ASSERT_NO_THROW(DAGMC::OBBUtils::constructBasisCont(elems, basis, points));
  ASSERT_EQ(points.n_cols, nNodesPerFace);
  ASSERT_EQ(points.n_rows, 3);
  orthonormal = checkBasis(basis, tol, errmsg);
  EXPECT_TRUE(orthonormal) << errmsg.str();

}

TEST_F(OBBTetTest, OBBConstructor) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);

  const DAGMC::ConstructMethod methods[]
    = { DAGMC::ConstructMethod::CONT,
        DAGMC::ConstructMethod::DISCRETE,
        DAGMC::ConstructMethod::ALIGNED
      };

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elEnd
    = meshPtr->elements_end();

  std::stringstream err;
  for (auto& method : methods) {
    err << "Failed for method " << method << "\n";
    DAGMC::OrientedBoundingBox obb(elBeg, elEnd, method);
    EXPECT_TRUE(obb.isConstructed()) << err.str();
    EXPECT_TRUE(obb.isSane()) << err.str();
    EXPECT_EQ(obb.status(), DAGMC::Box::success);
    EXPECT_NO_THROW(obb.getBox()) << err.str();
    for (auto& p : pointsLM) {
      EXPECT_TRUE(obb.containsPoint(p)) << err.str() << p;
    }
  }

}


TEST_F(OBBTetTest, SharedPtrTest) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);

  DAGMC::const_element_iterator elBeg
    = meshPtr->elements_begin();
  DAGMC::const_element_iterator elNext = elBeg;
  elNext++;

  // To use shared_from_this functionality, obb MUST
  // created as a shared ptr
  std::shared_ptr<DAGMC::OrientedBoundingBox> obbptr =
      std::make_shared<DAGMC::OrientedBoundingBox>(elBeg, elNext);

  EXPECT_EQ(obbptr.use_count(), 1);

  {
    // Test cast
    std::shared_ptr<DAGMC::TreeNode> node = obbptr->shared_from_this();
    EXPECT_NE(node, nullptr);
    EXPECT_EQ(node.use_count(), 2);
    EXPECT_EQ(obbptr.use_count(), 2);
  }

  EXPECT_EQ(obbptr.use_count(), 1);

  // Bad usage
  DAGMC::OrientedBoundingBox obb(elBeg, elNext);
  EXPECT_THROW(obb.shared_from_this(), std::bad_weak_ptr);

}

//---------------------------------------------------------------------------//

TEST_F(OBBFileTest, Read) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);

  for (unsigned int i = 0; i < files.size(); i++) {
    //Attempt to read file
    EXPECT_TRUE(Read(files.at(i)));
    EXPECT_TRUE(getSurfsAndElems());
    EXPECT_FALSE(libMeshException);

    // Move onto next file if we failed to read this one
    if (libMeshException)
      continue;

    //Check that mesh has the right number of surfs, elems, nodes.
    EXPECT_EQ(nSurfs(), 6);
    for (unsigned int iSurf = 0; iSurf < nSurfs(); iSurf++) {
      EXPECT_EQ(nElemsInSurf(iSurf), 26);
      EXPECT_EQ(nodesPerElemInSurf(iSurf), 3);
    }
  }
}

TEST_F(OBBFileTest, OBBUtils) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);

  // Loop over files to test
  std::stringstream err;
  for (unsigned int i = 0; i < files.size(); i++) {
    err.str("");
    err << "Test failed for file " << files.at(i) << "\n";
    // Don't continue if we fail read.
    ASSERT_TRUE(Read(files.at(i))) << err.str();
    ASSERT_FALSE(libMeshException) << err.str();

    DAGMC::const_element_iterator elBeg
      = meshPtr->elements_begin();
    DAGMC::const_element_iterator elEnd
      = meshPtr->elements_end();
    DAGMC::ElemConstItContainer elems(elBeg, elEnd);

    DAGMC::Matrix points;
    ASSERT_NO_THROW(DAGMC::OBBUtils::getPointsMatrix(elems, points)) << err.str();
    // 468 = 6 faces * 26 elems / face * 3 nodes /elem
    EXPECT_EQ(points.n_cols, 468) << err.str();
    EXPECT_EQ(points.n_rows, 3) << err.str();

    // Find the min/max points along x,y,z
    DAGMC::Matrix basis = { {1., 0., 0}, {0., 1., 0.}, {0., 0., 1.} };
    DAGMC::Vector min;
    DAGMC::Vector max;
    ASSERT_NO_THROW(DAGMC::OBBUtils::findExtremalPoints(points, basis, min, max)) << err.str();
    EXPECT_EQ(min.n_rows, 3) << err.str();
    EXPECT_EQ(max.n_rows, 3) << err.str();
    ASSERT_EQ(min.n_rows, max.n_rows) << err.str();

    for (unsigned int irow = 0; irow < max.n_rows; irow++) {
      std::stringstream ss;
      ss << "row " << irow;
      double diff = fabs(min(irow) + vals.at(i));
      EXPECT_LT(diff, tol) << err.str() << ss.str();
      diff = fabs(max(irow) - vals.at(i));
      EXPECT_LT(diff, tol) << err.str() << ss.str();;
    }

    // Find the basis from the points' covariance
    basis.reset();
    DAGMC::Matrix covMat = arma::cov(points.t());
    EXPECT_EQ(covMat.n_rows, 3) << err.str();
    EXPECT_EQ(covMat.n_cols, 3) << err.str();
    EXPECT_TRUE(covMat.is_symmetric()) << err.str();

    ASSERT_NO_THROW(DAGMC::OBBUtils::constructBasisFromCov(covMat, basis)) << err.str();
    std::stringstream msg;
    bool orthonormal = checkBasis(basis, tol, msg);
    EXPECT_TRUE(orthonormal) << err.str() << msg.str();

    // Find element statistics
    points.reset();
    std::vector<double> areas;
    DAGMC::Vector mean;
    ASSERT_NO_THROW(DAGMC::OBBUtils::getElemStats(elems, areas, mean, points)) << err.str();
    // 468 = 6 faces * 26 elems / face * 3 nodes /elem
    EXPECT_EQ(points.n_cols, 468) << err.str();
    // 156 = 6 faces * 26 elems / face
    EXPECT_EQ(areas.size(), 156) << err.str();
    EXPECT_EQ(points.n_rows, 3) << err.str();
    EXPECT_EQ(mean.n_rows, 3) << err.str();

    // Check mean makes sense
    for (unsigned int irow = 0; irow < max.n_rows; irow++) {
      EXPECT_LT(mean.at(irow), vals.at(i)) << err.str();
      EXPECT_GT(mean.at(irow), -vals.at(i)) << err.str();
    }

    // Compare areas to element areas
    DAGMC::const_element_iterator elIt = elBeg;
    unsigned iElem = 0;
    for (; elIt != elEnd; ++elIt) {
      ASSERT_LT(iElem, areas.size());
      double area = areas.at(iElem);
      double areaTest = (**elIt).volume();
      double diff = fabs(area - areaTest);
      EXPECT_LT(diff, tol);
      iElem++;
    }

    // Get covariance matrix from element stats
    covMat.reset();
    ASSERT_NO_THROW(DAGMC::OBBUtils::calcCov(areas, mean, points, covMat)) << err.str();
    EXPECT_EQ(covMat.n_rows, 3) << err.str();
    EXPECT_EQ(covMat.n_cols, 3) << err.str();
    EXPECT_TRUE(covMat.is_symmetric()) << err.str();

  }
}

TEST_F(OBBFileTest, OBBconstructor) {

  // Don't continue if libMesh threw an exception
  ASSERT_FALSE(libMeshException);
  // Don't continue if we have null pointers
  ASSERT_NE(meshPtr, nullptr);

  // Loop over files to test
  std::stringstream err;
  for (unsigned int i = 0; i < files.size(); i++) {
    err.str("");
    err << "Test failed for file " << files.at(i) << "\n";
    // Don't continue if we fail read.
    ASSERT_TRUE(Read(files.at(i))) << err.str();
    ASSERT_FALSE(libMeshException) << err.str();

    {
      DAGMC::const_element_iterator elBeg
        = meshPtr->elements_begin();
      DAGMC::const_element_iterator elEnd
        = meshPtr->elements_end();
      DAGMC::OrientedBoundingBox obb(elBeg, elEnd);

      bool constructed = obb.isConstructed();
      EXPECT_TRUE(constructed) << err.str();
      if (constructed) {
        EXPECT_TRUE(obb.isSane()) << err.str();
        EXPECT_EQ(obb.status(), DAGMC::Box::success) << err.str();

        // Check box contains all nodes in mesh
        unsigned int nNodes = meshPtr->n_nodes();
        for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
          const libMesh::Node* nodePtr =  meshPtr->node_ptr(0);
          ASSERT_NE(nodePtr, nullptr);
          EXPECT_TRUE(obb.containsPoint(*nodePtr));
        }
      }
    }

    // Loop over surfaces
    for (auto& isurf : surfIDs) {
      err.str("");
      err << "Test failed for file " << files.at(i)
          << " surface " << isurf << "\n";
      // Fetch all elements in mesh
      DAGMC::const_element_iterator elBeg
        = meshPtr->active_subdomain_elements_begin(isurf);;
      DAGMC::const_element_iterator elEnd
        = meshPtr->active_subdomain_elements_end(isurf);
      // Make an obb for this surface
      DAGMC::OrientedBoundingBox obb(elBeg, elEnd);
      EXPECT_TRUE(obb.isConstructed()) << err.str() ;
      EXPECT_TRUE(obb.isSane()) << err.str();

      // Expect degenerate along one direction
      ASSERT_NO_THROW(obb.getBox());
      const DAGMC::Box& box = obb.getBox();
      EXPECT_EQ(box.nDegenerate(), 1);
    }

  }
}

//---------------------------------------------------------------------------//

bool checkBasis(const DAGMC::Matrix& basis, double tol, std::stringstream& errmsg) {

  errmsg.str("");
  if (basis.n_rows != 3) {
    errmsg << "Basis failed row number test. N rows = ";
    errmsg << basis.n_rows;
    return false;
  }
  if (basis.n_cols != 3) {
    errmsg << "Basis failed col number test. N cols = ";
    errmsg << basis.n_cols;
    return false;
  }
  for (unsigned int icol = 0; icol < basis.n_cols; icol++) {
    // Check normalisation
    const DAGMC::Vector& ivec = basis.col(icol);
    double vnorm = arma::norm(ivec);
    if (fabs(vnorm - 1.0) > tol) {
      errmsg << "Basis failed normalisation test. Norm vec ";
      errmsg << icol << " = " << vnorm;
      return false;
    }
    for (unsigned int jcol = icol + 1; jcol < basis.n_cols; jcol++) {
      //Check orthogonal
      const DAGMC::Vector& jvec = basis.col(jcol);
      double prod = arma::dot(ivec, jvec);
      if (fabs(vnorm - 1.0) > tol) {
        errmsg << "Basis failed orthogonality between test between vecs ";
        errmsg << icol << "," << jcol;
        return false;
      }
    }
  }
  return true;
};
