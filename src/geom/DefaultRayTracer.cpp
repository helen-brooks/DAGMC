#include "DefaultRayTracer.hpp"

namespace DAGMC {

DefaultRayTracer::DefaultRayTracer(std::shared_ptr<MoabInterface> mesh_interface) {

  // Get the geometry topo tool
  GTT  = mesh_interface->gtt();

#ifdef DOUBLE_DOWN
  std::cout << "Using the DOUBLE-DOWN interface to Embree." << std::endl;
  RTI = std::make_unique<RayTracingInterface>(GTT);
#else
  RTI = std::make_unique<RayTracingInterface>(GTT.get());
#endif

  errHandler = std::make_unique<MoabErrHandler>();

}

ErrorCode DefaultRayTracer::init() {

  // find all geometry sets
  errHandler->checkSetErr(GTT->find_geomsets(),
                          "Could not find the geometry sets");

  // implicit complement
  errHandler->checkSetErr(init_implicit_complement(),
                          "Failed to setup the implicit complement");

  // build obbs
  errHandler->checkSetErr(init_obb(), "Failed to setup the OBBs");

  return DAG_SUCCESS;
}

ErrorCode DefaultRayTracer::init_implicit_complement() {
  // If it doesn't already exist, create implicit complement
  // Create data structures for implicit complement
  ErrorCode rval = ErrorCode(GTT->setup_implicit_complement());
  if (DAG_SUCCESS != rval) {
    std::cerr << "Failed to find or create implicit complement handle." << std::endl;
    return rval;
  }
  return DAG_SUCCESS;
}

ErrorCode DefaultRayTracer::init_obb()

{
  // If we havent got an OBB Tree, build one.
  if (!GTT->have_obb_tree()) {
    std::cout << "Building acceleration data structures..." << std::endl;
#ifdef DOUBLE_DOWN
    errHandler->checkSetErr(RTI->init(),
                            "Failed to build obb trees");
#else
    errHandler->checkSetErr(GTT->construct_obb_trees(),
                            "Failed to build obb trees");
#endif
  }
  return DAG_SUCCESS;
}

ErrorCode DefaultRayTracer::next_vol(EntityHandle surface, EntityHandle old_volume,
                                     EntityHandle& new_volume) {
  return ErrorCode(GTT->next_vol(surface, old_volume, new_volume));
}

bool DefaultRayTracer::is_implicit_complement(EntityHandle volume) {
  return ErrorCode(GTT->is_implicit_complement(volume));
}

// get sense of surface(s) wrt volume
ErrorCode DefaultRayTracer::surface_sense(EntityHandle volume, int num_surfaces,
                                          const EntityHandle* surfaces, int* senses_out) {
  return ErrorCode(GTT->get_surface_senses(volume, num_surfaces, surfaces,
                                           senses_out));
}

// get sense of surface(s) wrt volume
ErrorCode DefaultRayTracer::surface_sense(EntityHandle volume, EntityHandle surface,
                                          int& sense_out) {
  return ErrorCode(GTT->get_sense(surface, volume, sense_out));
}

ErrorCode DefaultRayTracer::point_in_volume(const moab::EntityHandle volume,
                                            const double xyz[3],
                                            int& result,
                                            const double* uvw,
                                            const moab::GeomQueryTool::RayHistory* history) {
  return ErrorCode(RTI->point_in_volume(volume, xyz, result, uvw, history));
}

ErrorCode DefaultRayTracer::point_in_volume_slow(moab::EntityHandle volume,
                                                 const double xyz[3],
                                                 int& result) {
  return ErrorCode(RTI->point_in_volume_slow(volume, xyz, result));
}

ErrorCode DefaultRayTracer::test_volume_boundary(const moab::EntityHandle volume,
                                                 const moab::EntityHandle surface,
                                                 const double xyz[3],
                                                 const double uvw[3],
                                                 int& result,
                                                 const moab::GeomQueryTool::RayHistory* history) {
  return ErrorCode(RTI->test_volume_boundary(volume, surface, xyz, uvw, result, history));
}

ErrorCode DefaultRayTracer::ray_fire(const moab::EntityHandle volume,
                                     const double point[3],
                                     const double dir[3],
                                     moab::EntityHandle& next_surf,
                                     double& next_surf_dist,
                                     moab::GeomQueryTool::RayHistory* history,
                                     double user_dist_limit,
                                     int ray_orientation,
                                     OrientedBoxTreeTool::TrvStats* stats) {
  return ErrorCode(RTI->ray_fire(volume, point, dir, next_surf, next_surf_dist,
                                 history, user_dist_limit, ray_orientation,
                                 stats));
}

ErrorCode DefaultRayTracer::closest_to_location(moab::EntityHandle volume,
                                                const double point[3],
                                                double& result,
                                                moab::EntityHandle* closest_surf) {
  return ErrorCode(RTI->closest_to_location(volume, point, result, closest_surf));
}

ErrorCode DefaultRayTracer::get_normal(moab::EntityHandle surf,
                                       const double loc[3],
                                       double angle[3],
                                       const moab::GeomQueryTool::RayHistory* history) {
  return ErrorCode(RTI->get_normal(surf, loc, angle, history));
}

ErrorCode DefaultRayTracer::measure_volume(moab::EntityHandle volume,
                                           double& result) {
  return ErrorCode(RTI->measure_volume(volume, result));
}

ErrorCode DefaultRayTracer::measure_area(moab::EntityHandle surface,
                                         double& result) {
  return ErrorCode(RTI->measure_area(surface, result));
}

double DefaultRayTracer::get_numerical_precision() {
  return RTI->get_numerical_precision();
}

double DefaultRayTracer::get_overlap_thickness() {
  return RTI->get_overlap_thickness();
}

void DefaultRayTracer::set_numerical_precision(double val) {
  RTI->set_numerical_precision(val);
}

void DefaultRayTracer::set_overlap_thickness(double val) {
  RTI->set_overlap_thickness(val);
}


}
