#include "DagMCmoab.hpp"

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

DagMCmoab::DagMCmoab(std::shared_ptr<moab::Interface> mb_impl, double overlap_tolerance, double p_numerical_precision) {

  // Create an interface to MOAB
  moab_interface = std::make_shared<MoabInterface>(mb_impl);

  init(overlap_tolerance, p_numerical_precision);
}

DagMCmoab::DagMCmoab(moab::Interface* mb_impl, double overlap_tolerance, double p_numerical_precision) {

  // Create an interface to MOAB
  moab_interface = std::make_shared<MoabInterface>(mb_impl);

  init(overlap_tolerance, p_numerical_precision);

}

void DagMCmoab::init(double overlap_tolerance, double numerical_precision) {

  // Copy moab interface pointer into generic base class pointer
  mesh_interface = moab_interface;

  // Create error handler
  errHandler = std::make_unique<MoabErrHandler>();

  ray_tracer = std::make_unique<DefaultRayTracer>(moab_interface, overlap_tolerance, numerical_precision);

}

// *****************************************************************************
// SECTION III
// *****************************************************************************

EntityHandle DagMCmoab::entity_by_id(int dimension, int id) {
  return moab_interface->entity_by_id(dimension, id);
}

int DagMCmoab::get_entity_id(EntityHandle this_ent) {
  return moab_interface->get_entity_id(this_ent);
}

EntityHandle DagMCmoab::entity_by_index(int dimension, int index) {
  return moab_interface->entity_by_index(dimension, index);
}

int DagMCmoab::index_by_handle(EntityHandle handle) {
  return moab_interface->index_by_handle(handle);
}

int DagMCmoab::id_by_index(int dimension, int index) {
  return moab_interface->id_by_index(dimension, index);
}

unsigned int DagMCmoab::num_entities(int dimension) {
  return moab_interface->num_entities(dimension);
}

// *****************************************************************************
// SECTION IV: Handling DagMC settings
// *****************************************************************************

double DagMCmoab::faceting_tolerance()  {
  return moab_interface->get_faceting_tol();
};

// *****************************************************************************
// SECTION V: Metadata handling
// *****************************************************************************

ErrorCode DagMCmoab::detect_available_props(std::vector<std::string>& keywords_list,
                                            const char* delimiters) {

  if (!moab_interface->get_keywords(keywords_list, delimiters)) {
    return moab_interface->code();
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
  if (!moab_interface->update_properties(prop_names, delimiters)) {
    return moab_interface->code();
  }

  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::prop_value(EntityHandle eh, const std::string& prop, std::string& value) {

  if (!moab_interface->get_property(eh, prop, value)) {
    return moab_interface->code();
  }

  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::prop_values(EntityHandle eh, const std::string& prop,
                                 std::vector< std::string >& values) {

  if (!moab_interface->get_properties(eh, prop, values)) {
    return moab_interface->code();
  }

  return DAG_SUCCESS;

}

bool DagMCmoab::has_prop(EntityHandle eh, const std::string& prop) {

  return moab_interface->has_property(eh, prop);
}

// TO-DO: figure out where this are used
ErrorCode DagMCmoab::get_all_prop_values(const std::string& prop, std::vector<std::string>& return_list) {

  std::set<EntityHandle> dummy;
  std::set<std::string> unique_vals;
  if (!moab_interface->get_ents_and_vals_with_prop(prop, dummy, unique_vals)) {
    return moab_interface->code();
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
  if (!moab_interface->get_ents_and_vals_with_prop(prop, handles, dummy, checkval, dimension, strval)) {
    return moab_interface->code();
  }
  return_list.assign(handles.begin(), handles.end());
  return DAG_SUCCESS;
}

} // namespace DAGMC
