#include "DagMCmoab.hpp"

#ifdef DOUBLE_DOWN
#include "RTI.hpp"
#include "MOABRay.h"
#endif

//#define MB_OBB_TREE_TAG_NAME "OBB_TREE"

namespace DAGMC {

/* Tolerance Summary

   Facet Tolerance:
   Maximum distance between continuous solid model surface and faceted surface.
     Performance:  increasing tolerance increased performance (fewer triangles)
     Robustness:   should not be affected
     Knowledge:    user must understand how coarser faceting influences accuracy
                   of results
*/

const bool counting = false; /* controls counts of ray casts and pt_in_vols */

DagMCmoab::DagMCmoab(std::shared_ptr<Interface> mb_impl, double overlap_tolerance, double p_numerical_precision) {

#ifdef DOUBLE_DOWN
  std::cout << "Using the DOUBLE-DOWN interface to Embree." << std::endl;
#endif

  // Create error handler
  errHandler = std::make_unique<MoabErrHandler>();

  // Create an interface to MOAB
  mesh_interface = std::make_shared<MoabInterface>(mb_impl);

  // Get the geometry topo tool
  GTT  = mesh_interface->gtt();

#ifdef DOUBLE_DOWN
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT));
#else
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT.get()));
#endif
  this->set_overlap_thickness(overlap_tolerance);
  this->set_numerical_precision(p_numerical_precision);

}

DagMCmoab::DagMCmoab(Interface* mb_impl, double overlap_tolerance, double p_numerical_precision) {

  // Create error handler
  errHandler = std::make_unique<MoabErrHandler>();

  // Create an interface to MOAB
  mesh_interface = std::make_shared<MoabInterface>(mb_impl);

  // Get the geometry topo tool
  GTT  = mesh_interface->gtt();

#ifdef DOUBLE_DOWN
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT));
#else
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT.get()));
#endif
  this->set_overlap_thickness(overlap_tolerance);
  this->set_numerical_precision(p_numerical_precision);

}

// *****************************************************************************
// SECTION I: Geometry Initialization and problem setup
// *****************************************************************************

// the standard DAGMC load file method
ErrorCode DagMCmoab::load_file(const char* cfile) {

  if (!mesh_interface->load(std::string(cfile))) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;

}

// helper function to load the existing contents of a MOAB instance into DAGMC
ErrorCode DagMCmoab::load_existing_contents() {

  if (!mesh_interface->setup_geom()) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;

}

// Set up the implicit complement
ErrorCode DagMCmoab::setup_impl_compl() {
  // If it doesn't already exist, create implicit complement
  // Create data structures for implicit complement
  ErrorCode rval = ErrorCode(GTT->setup_implicit_complement());
  if (DAG_SUCCESS != rval) {
    std::cerr << "Failed to find or create implicit complement handle." << std::endl;
    return rval;
  }
  return DAG_SUCCESS;
}

// initialise the obb tree
ErrorCode DagMCmoab::init_OBBTree() {

  // find all geometry sets
  errHandler->checkSetErr(GTT->find_geomsets(),
                          "Could not find the geometry sets");

  // implicit complement
  errHandler->checkSetErr(setup_impl_compl(),
                          "Failed to setup the implicit complement");

  // build obbs
  errHandler->checkSetErr(setup_obbs(), "Failed to setup the OBBs");

  // setup indices
  errHandler->checkSetErr(setup_indices(), "Failed to setup problem indices");

  return DAG_SUCCESS;
}

// setups of the indices for the problem, builds a list of surface and volumes
// indices
ErrorCode DagMCmoab::setup_indices() {

  if (!mesh_interface->setup_indices()) {
    errHandler->checkSetErr(mesh_interface->code(),
                            "Failed to build surface/volume indices");
  }
  return DAG_SUCCESS;
}

// sets up the obb tree for the problem
ErrorCode DagMCmoab::setup_obbs() {
  // If we havent got an OBB Tree, build one.
  if (!GTT->have_obb_tree()) {
    std::cout << "Building acceleration data structures..." << std::endl;
#ifdef DOUBLE_DOWN
    errHandler->checkSetErr(ray_tracer->init(),
                            "Failed to build obb trees");
#else
    errHandler->checkSetErr(GTT->construct_obb_trees(),
                            "Failed to build obb trees");
#endif
  }
  return DAG_SUCCESS;
}

// *****************************************************************************
// SECTION II: Fundamental Geometry Operations/Queries
// *****************************************************************************

ErrorCode DagMCmoab::ray_fire(const EntityHandle volume, const double point[3],
                              const double dir[3], EntityHandle& next_surf,
                              double& next_surf_dist,
                              RayHistory* history,
                              double user_dist_limit, int ray_orientation,
                              OrientedBoxTreeTool::TrvStats* stats) {
  return ErrorCode(ray_tracer->ray_fire(volume, point, dir, next_surf, next_surf_dist,
                                        history, user_dist_limit, ray_orientation,
                                        stats));
}

ErrorCode DagMCmoab::point_in_volume(const EntityHandle volume, const double xyz[3],
                                     int& result, const double* uvw,
                                     const RayHistory* history) {
  return ErrorCode(ray_tracer->point_in_volume(volume, xyz, result, uvw, history));
}

ErrorCode DagMCmoab::test_volume_boundary(const EntityHandle volume,
                                          const EntityHandle surface,
                                          const double xyz[3], const double uvw[3],
                                          int& result,
                                          const RayHistory* history) {
  return ErrorCode(ray_tracer->test_volume_boundary(volume, surface, xyz, uvw, result,
                                                    history));
}

// use spherical area test to determine inside/outside of a polyhedron.
ErrorCode DagMCmoab::point_in_volume_slow(EntityHandle volume, const double xyz[3],
                                          int& result) {
  return ErrorCode(ray_tracer->point_in_volume_slow(volume, xyz, result));
}

// detemine distance to nearest surface
ErrorCode DagMCmoab::closest_to_location(EntityHandle volume,
                                         const double coords[3], double& result,
                                         EntityHandle* surface) {
  return ErrorCode(ray_tracer->closest_to_location(volume, coords, result, surface));
}

// calculate volume of polyhedron
ErrorCode DagMCmoab::measure_volume(EntityHandle volume, double& result) {
  return ErrorCode(ray_tracer->measure_volume(volume, result));
}

// sum area of elements in surface
ErrorCode DagMCmoab::measure_area(EntityHandle surface, double& result) {
  return ErrorCode(ray_tracer->measure_area(surface, result));

}

// get sense of surface(s) wrt volume
ErrorCode DagMCmoab::surface_sense(EntityHandle volume, int num_surfaces,
                                   const EntityHandle* surfaces, int* senses_out) {
  return ErrorCode(GTT->get_surface_senses(volume, num_surfaces, surfaces,
                                           senses_out));
}

// get sense of surface(s) wrt volume
ErrorCode DagMCmoab::surface_sense(EntityHandle volume, EntityHandle surface,
                                   int& sense_out) {
  return ErrorCode(GTT->get_sense(surface, volume, sense_out));
}

ErrorCode DagMCmoab::get_angle(EntityHandle surf, const double in_pt[3],
                               double angle[3],
                               const RayHistory* history) {
  return ErrorCode(ray_tracer->get_normal(surf, in_pt, angle, history));
}

ErrorCode DagMCmoab::next_vol(EntityHandle surface, EntityHandle old_volume,
                              EntityHandle& new_volume) {
  return ErrorCode(GTT->next_vol(surface, old_volume, new_volume));
}

bool DagMCmoab::is_implicit_complement(EntityHandle volume) {
  return GTT->is_implicit_complement(volume);
}

// *****************************************************************************
// SECTION III
// *****************************************************************************

EntityHandle DagMCmoab::entity_by_id(int dimension, int id) {
  return mesh_interface->entity_by_id(dimension, id);
}

int DagMCmoab::get_entity_id(EntityHandle this_ent) {
  return mesh_interface->get_entity_id(this_ent);
}

EntityHandle DagMCmoab::entity_by_index(int dimension, int index) {
  return mesh_interface->entity_by_index(dimension, index);
}

int DagMCmoab::index_by_handle(EntityHandle handle) {
  return mesh_interface->index_by_handle(handle);
}

int DagMCmoab::id_by_index(int dimension, int index) {
  return mesh_interface->id_by_index(dimension, index);
}

unsigned int DagMCmoab::num_entities(int dimension) {
  return mesh_interface->num_entities(dimension);
}

// *****************************************************************************
// SECTION IV
// *****************************************************************************

double DagMCmoab::overlap_thickness() {
  return ray_tracer->get_overlap_thickness();
}

double DagMCmoab::numerical_precision() {
  return ray_tracer->get_numerical_precision();
}

void DagMCmoab::set_overlap_thickness(double new_thickness) {
  ray_tracer->set_overlap_thickness(new_thickness);
}

void DagMCmoab::set_numerical_precision(double new_precision) {
  ray_tracer->set_numerical_precision(new_precision);

}

// *****************************************************************************
// SECTION V: Metadata handling
// *****************************************************************************

ErrorCode DagMCmoab::detect_available_props(std::vector<std::string>& keywords_list,
                                            const char* delimiters) {

  if (!mesh_interface->get_keywords(keywords_list, delimiters)) {
    return mesh_interface->code();
  }
  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::parse_properties(const std::vector<std::string>& keywords,
                                      const std::map<std::string,
                                      std::string>& keyword_synonyms,
                                      const char* delimiters) {

  // Master keyword map, mapping user-set words in cubit to canonical property names
  std::map< std::string, std::string > keyword_map(keyword_synonyms);

  // Now add the requested keywords
  for (auto key : keywords) {
    keyword_map[key] = key;
  }

  // Create the set of all canonical property names
  std::set< std::string > prop_names;
  for (auto keypair : keyword_map) {
    prop_names.insert(keypair.second);
  }

  // Set up DagMC's property tags based on what's been requested
  if (!mesh_interface->update_properties(prop_names, delimiters)) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::prop_value(EntityHandle eh, const std::string& prop, std::string& value) {

  if (!mesh_interface->get_property(eh, prop, value)) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::prop_values(EntityHandle eh, const std::string& prop,
                                 std::vector< std::string >& values) {

  if (!mesh_interface->get_properties(eh, prop, values)) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;

}

bool DagMCmoab::has_prop(EntityHandle eh, const std::string& prop) {

  return mesh_interface->has_property(eh, prop);
}

// TO-DO: figure out where this are used
ErrorCode DagMCmoab::get_all_prop_values(const std::string& prop, std::vector<std::string>& return_list) {

  std::set<EntityHandle> dummy;
  std::set<std::string> unique_vals;
  if (!mesh_interface->get_ents_and_vals_with_prop(prop, dummy, unique_vals)) {
    return mesh_interface->code();
  }
  return_list.assign(unique_vals.begin(), unique_vals.end());
  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::entities_by_property(const std::string& prop,
                                          std::vector<EntityHandle>& return_list,
                                          int dimension, const std::string* value) {

  std::string strval = "";
  bool checkval = false;
  if (value != nullptr) {
    checkval = true;
    strval = *value;
  }

  std::set<EntityHandle> handles;
  std::set<std::string> dummy;
  if (!mesh_interface->get_ents_and_vals_with_prop(prop, handles, dummy, checkval, dimension, strval)) {
    return mesh_interface->code();
  }
  return_list.assign(handles.begin(), handles.end());
  return DAG_SUCCESS;
}

// ***************************************************************************
// SECTION VI
// ***************************************************************************

ErrorCode DagMCmoab::write_mesh(const char* ffile) {

  // Write out a mesh file if requested
  if (ffile) {
    if (!mesh_interface->write(std::string(ffile))) {
      std::cerr << "Failed to write mesh to " << ffile << "." << std::endl;
      return mesh_interface->code();
    }
  }

  return DAG_SUCCESS;
}


} // namespace DAGMC
