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
  errHandler->checkSetErr(find_geomsets(),
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
  return GTT->entity_by_id(dimension, id);
}

int DagMCmoab::get_entity_id(EntityHandle this_ent) {
  return GTT->global_id(this_ent);
}

// *****************************************************************************
// SECTION IV
// *****************************************************************************

double DagMCmoab::overlap_thickness() { return ray_tracer->get_overlap_thickness(); }

double DagMCmoab::numerical_precision() { return ray_tracer->get_numerical_precision(); }

void DagMCmoab::set_overlap_thickness(double new_thickness) {
  ray_tracer->set_overlap_thickness(new_thickness);
}

void DagMCmoab::set_numerical_precision(double new_precision) {
  ray_tracer->set_numerical_precision(new_precision);

}

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

  for (std::vector<std::string>::const_iterator i = keywords.begin();
       i != keywords.end(); ++i) {
    keyword_map[*i] = *i;
  }

  // The set of all canonical property names
  std::set< std::string > prop_names;
  for (auto i = keyword_map.begin();
       i != keyword_map.end(); ++i) {
    prop_names.insert((*i).second);
  }

  // // Now add the requested keywords
  // for (auto key : keywords){
  //   keyword_map[key] = key;
  // }

  // // Create the set of all canonical property names
  // std::set< std::string > prop_names;
  // for (auto keypair : keyword_map) {
  //   prop_names.insert(keypair.second);
  // }

  // Set up DagMC's property tags based on what's been requested
  if (!mesh_interface->set_tagmap(prop_names, property_tagmap)) {
    return mesh_interface->code();
  }

  // Now that tags are ready, append moab's group properties
  if (!mesh_interface->append_group_properties(delimiters, property_tagmap)) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::prop_value(EntityHandle eh, const std::string& prop, std::string& value) {

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return DAG_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;

  if (!mesh_interface->get_tag_name(proptag, eh, value)) {
    return mesh_interface->code();
  }

  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::prop_values(EntityHandle eh, const std::string& prop,
                                 std::vector< std::string >& values) {

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return DAG_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;

  // Convert a property tag's value on a handle to a list of strings
  if (!mesh_interface->get_tag_data_vec(proptag, eh, values)) {
    return mesh_interface->code();
  }
  return DAG_SUCCESS;

}

bool DagMCmoab::has_prop(EntityHandle eh, const std::string& prop) {

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return false;
  }
  Tag proptag = (*it).second;
  std::string dummyvalue;
  bool found = mesh_interface->get_tag_name(proptag, eh, dummyvalue);

  return found;

}

// TO-DO: figure out where this are used
ErrorCode DagMCmoab::get_all_prop_values(const std::string& prop, std::vector<std::string>& return_list) {

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return DAG_TAG_NOT_FOUND;
  }

  EntityHandle root = 0;
  Tag proptag = (*it).second;
  Range all_ents;
  if (!mesh_interface->get_tagged_entity_sets(root, proptag, all_ents)) {
    return mesh_interface->code();
  }

  std::set<std::string> unique_values;
  for (auto& entity : all_ents) {
    std::vector<std::string> values;
    ErrorCode rval = rval = prop_values(entity, prop, values);
    if (DAG_SUCCESS != rval)
      return rval;
    unique_values.insert(values.begin(), values.end());
  }

  return_list.assign(unique_values.begin(), unique_values.end());
  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::entities_by_property(const std::string& prop,
                                          std::vector<EntityHandle>& return_list,
                                          int dimension, const std::string* value) {
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return DAG_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;
  EntityHandle root = 0;
  std::vector<Tag> tags = {proptag, GTT->get_geom_tag()};

  // Note that we cannot specify values for proptag here-- the passed value,
  // if it exists, may be only a subset of the packed string representation
  // of this tag.

  if (!mesh_interface->get_tagged_entity_sets(root, tags, all_ents)) {
    return mesh_interface->code();
  }

  std::set<EntityHandle> handles;
  for (auto& ent : all_ents) {
    std::vector<std::string> values;
    rval = prop_values(ent, prop, values);
    if (DAG_SUCCESS != rval)
      return rval;
    if (value && std::find(values.begin(), values.end(), *value) != values.end()) {
      handles.insert(ent);
    } else {
      handles.insert(ent);
    }
  }

  return_list.assign(handles.begin(), handles.end());
  return DAG_SUCCESS;
}


} // namespace DAGMC
