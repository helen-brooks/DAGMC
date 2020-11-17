// This file was created by H Brooks on 29/06/2020
#include "DagMCBase.hpp"

using namespace DAGMC;

// *****************************************************************************
// Public methods
// *****************************************************************************

// get the float verision of dagmc version string
float DagMCBase::version(std::string* version_string) {
  if (NULL != version_string)
    *version_string = std::string("DagMC version ") + std::string(DAGMC_VERSION_STRING);
  return DAGMC_VERSION;
}


// *****************************************************************************
// SECTION I: Geometry Initialization and problem setup
// *****************************************************************************

// the standard DAGMC load file method
ErrorCode DagMCBase::load_file(const char* cfile) {

  if (!mesh_interface->load(std::string(cfile))) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;

}

// Helper function to load the existing contents of a already open file
ErrorCode DagMCBase::load_existing_contents() {

  if (!mesh_interface->finish_setup()) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;

}

// initialise the obb tree
ErrorCode DagMCBase::init_OBBTree() {

  // Find all geometry sets
  // TO-DO: check if this call is every really needed
  if (!mesh_interface->setup_geom()) {
    errHandler->checkSetErr(mesh_interface->code(),
                            "Could not find the geometry sets");
  }

  // Setup ray tracer
  errHandler->checkSetErr(ray_tracer->init(), "Failed to initialise ray tracer");

  // setup indices
  errHandler->checkSetErr(setup_indices(), "Failed to setup problem indices");

  return DAG_SUCCESS;
}

// Set up the implicit complement
ErrorCode DagMCBase::setup_impl_compl() {
  return ray_tracer->init_implicit_complement();
}


// sets up the obb tree for the problem
ErrorCode DagMCBase::setup_obbs() {
  return ray_tracer->init_obb();
}

// setups of the indices for the problem, builds a list of surface and volumes
// indices
ErrorCode DagMCBase::setup_indices() {

  if (!mesh_interface->setup_indices()) {
    errHandler->checkSetErr(mesh_interface->code(),
                            "Failed to build surface/volume indices");
  }
  return DAG_SUCCESS;
}

// *****************************************************************************
// SECTION II: Fundamental Geometry Operations/Queries
// *****************************************************************************

ErrorCode DagMCBase::ray_fire(const EntityHandle volume, const double point[3],
                              const double dir[3], EntityHandle& next_surf,
                              double& next_surf_dist,
                              moab::GeomQueryTool::RayHistory* history,
                              double user_dist_limit, int ray_orientation,
                              moab::OrientedBoxTreeTool::TrvStats* stats) {
  return ray_tracer->ray_fire(volume, point, dir, next_surf, next_surf_dist,
                              history, user_dist_limit, ray_orientation,
                              stats);
}

ErrorCode DagMCBase::point_in_volume(const EntityHandle volume, const double xyz[3],
                                     int& result, const double* uvw,
                                     const moab::GeomQueryTool::RayHistory* history) {
  return ray_tracer->point_in_volume(volume, xyz, result, uvw, history);
}

ErrorCode DagMCBase::test_volume_boundary(const EntityHandle volume,
                                          const EntityHandle surface,
                                          const double xyz[3], const double uvw[3],
                                          int& result,
                                          const moab::GeomQueryTool::RayHistory* history) {
  return ray_tracer->test_volume_boundary(volume, surface, xyz, uvw, result,
                                          history);
}

// use spherical area test to determine inside/outside of a polyhedron.
ErrorCode DagMCBase::point_in_volume_slow(EntityHandle volume, const double xyz[3],
                                          int& result) {
  return ray_tracer->point_in_volume_slow(volume, xyz, result);
}

// detemine distance to nearest surface
ErrorCode DagMCBase::closest_to_location(EntityHandle volume,
                                         const double coords[3], double& result,
                                         EntityHandle* surface) {
  return ray_tracer->closest_to_location(volume, coords, result, surface);
}

// calculate volume of polyhedron
ErrorCode DagMCBase::measure_volume(EntityHandle volume, double& result) {
  return ray_tracer->measure_volume(volume, result);
}

// sum area of elements in surface
ErrorCode DagMCBase::measure_area(EntityHandle surface, double& result) {
  return ray_tracer->measure_area(surface, result);
}

// get sense of surface(s) wrt volume
ErrorCode DagMCBase::surface_sense(EntityHandle volume, int num_surfaces,
                                   const EntityHandle* surfaces, int* senses_out) {
  return ray_tracer->surface_sense(volume, num_surfaces, surfaces,
                                   senses_out);
}

// get sense of surface(s) wrt volume
ErrorCode DagMCBase::surface_sense(EntityHandle volume, EntityHandle surface,
                                   int& sense_out) {
  return ray_tracer->surface_sense(surface, volume, sense_out);
}

ErrorCode DagMCBase::get_angle(EntityHandle surf, const double in_pt[3],
                               double angle[3],
                               const moab::GeomQueryTool::RayHistory* history) {
  return ray_tracer->get_normal(surf, in_pt, angle, history);
}

ErrorCode DagMCBase::next_vol(EntityHandle surface, EntityHandle old_volume,
                              EntityHandle& new_volume) {
  return ray_tracer->next_vol(surface, old_volume, new_volume);
}

bool DagMCBase::is_implicit_complement(EntityHandle volume) {
  return ray_tracer->is_implicit_complement(volume);
}


// *****************************************************************************
// SECTION IV: Handling DagMC settings
// *****************************************************************************

double DagMCBase::overlap_thickness() {
  return ray_tracer->get_overlap_thickness();
}

double DagMCBase::numerical_precision() {
  return ray_tracer->get_numerical_precision();
}

void DagMCBase::set_overlap_thickness(double new_thickness) {
  ray_tracer->set_overlap_thickness(new_thickness);
}

void DagMCBase::set_numerical_precision(double new_precision) {
  ray_tracer->set_numerical_precision(new_precision);
}

// ***************************************************************************
// SECTION VI: misc
// ***************************************************************************

ErrorCode DagMCBase::write_mesh(const char* ffile) {

  // Write out a mesh file if requested
  if (ffile) {
    if (!mesh_interface->write(std::string(ffile))) {
      std::cerr << "Failed to write mesh to " << ffile << "." << std::endl;
      return mesh_interface->code();
    }
  }

  return DAG_SUCCESS;
}
