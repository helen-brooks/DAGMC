#include "DagMCmoab.hpp"

/* #include <algorithm>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h> */

#ifndef M_PI  /* windows */
# define M_PI 3.14159265358979323846
#endif

#ifdef DOUBLE_DOWN
#include "RTI.hpp"
#include "MOABRay.h"
#endif

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define FACETING_TOL_TAG_NAME "FACETING_TOL"
static const int null_delimiter_length = 1;

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
														  
  moab_instance_created = false;

  // Create error handler
  errHandler = std::make_unique<MoabErrHandler>();

  // If we aren't handed a moab instance create one
  if (nullptr == mb_impl) {
    mb_impl = std::make_shared<Core>();
    moab_instance_created = true;
  }

  MBI_shared_ptr = mb_impl;
  // set the internal moab pointer
  MBI = MBI_shared_ptr.get();

  // make new GeomTopoTool and GeomQueryTool
  GTT = std::make_shared<GeomTopoTool> (MBI, false);
#ifdef DOUBLE_DOWN
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT));
#else
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT.get()));
#endif
  this->set_overlap_thickness(overlap_tolerance);
  this->set_numerical_precision(p_numerical_precision);

}

DagMCmoab::DagMCmoab(Interface* mb_impl, double overlap_tolerance, double p_numerical_precision) {
  moab_instance_created = false;

  // Create error handler
  errHandler = std::make_unique<MoabErrHandler>();

  // set the internal moab pointer
  MBI = mb_impl;
  MBI_shared_ptr = nullptr;

  // make new GeomTopoTool and GeomQueryTool
  GTT = std::make_shared<GeomTopoTool> (MBI, false);
#ifdef DOUBLE_DOWN
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT));
#else
  ray_tracer = std::unique_ptr<RayTracer>(new RayTracer(GTT.get()));
#endif
  this->set_overlap_thickness(overlap_tolerance);
  this->set_numerical_precision(p_numerical_precision);

}

// Destructor
DagMCmoab::~DagMCmoab() {
  // if we created the moab instance
  // clear it
  if (moab_instance_created) {
    MBI->delete_mesh();
  }
}

// get the float verision of dagmc version string
float DagMC::version(std::string* version_string) {
  if (NULL != version_string)
    *version_string = std::string("DagMC version ") + std::string(DAGMC_VERSION_STRING);
  return DAGMC_VERSION;
}

// *****************************************************************************
// SECTION I: Geometry Initialization and problem setup
// *****************************************************************************

// the standard DAGMC load file method
ErrorCode DagMCmoab::load_file(const char* cfile) {
  ErrorCode rval;
  std::string filename(cfile);
  std::cout << "Loading file " << cfile << std::endl;
  // load options
  char options[120] = {0};
  std::string file_ext = "" ; // file extension

  // get the last 4 chars of file .i.e .h5m .sat etc
  int file_extension_size = 4;
  if (filename.size() > file_extension_size) {
    file_ext = filename.substr(filename.size() - file_extension_size);
  }
  EntityHandle file_set;
  rval = ErrorCode(MBI->create_meshset(moab::MESHSET_SET, file_set));
  if (DAG_SUCCESS != rval)
    return rval;

  rval = ErrorCode(MBI->load_file(cfile, &file_set, options, NULL, 0, 0));

  if (DAG_UNHANDLED_OPTION == rval) {
    // Some options were unhandled; this is common for loading h5m files.
    // Print a warning if an option was unhandled for a file that does not end in '.h5m'
    std::string filename(cfile);
    if (file_ext != ".h5m") {
      std::cerr << "DagMC warning: unhandled file loading options." << std::endl;
    }
  } else if (DAG_SUCCESS != rval) {
    std::cerr << "DagMC Couldn't read file " << cfile << std::endl;
    std::string message;
    if (DAG_SUCCESS == ErrorCode(MBI->get_last_error(message)) && !message.empty())
      std::cerr << "Error message: " << message << std::endl;

    return rval;
  }

  return finish_loading();
}

// setup the implicit compliment
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

// gets the entity sets tagged with geomtag 2 and 3
// surfaces and volumes respectively
ErrorCode DagMCmoab::setup_geometry(Range& surfs, Range& vols) {

  // get all surfaces
  errHandler->checkSetErr(GTT->get_gsets_by_dimension(2, surfs),
                          "Could not get surfaces from GTT");

  // get all volumes
  errHandler->checkSetErr(GTT->get_gsets_by_dimension(3, vols),
                          "Could not get volumes from GTT");

  return DAG_SUCCESS;
}

// sets up the obb tree for the problem
ErrorCode DagMCmoab::setup_obbs() {
  // If we havent got an OBB Tree, build one.
  if (!GTT->have_obb_tree()) {
    std::cout << "Building acceleration data structures..." << std::endl;
#ifdef DOUBLE_DOWN
    rval = 
    errHandler->checkSetErr(ray_tracer->init();,
                            "Failed to build obb trees");
#else
    errHandler->checkSetErr(GTT->construct_obb_trees(),
                            "Failed to build obb trees");
#endif
  }
  return DAG_SUCCESS;
}

// helper function to finish setting up required tags.
ErrorCode DagMCmoab::finish_loading() {
  ErrorCode rval;

  nameTag = get_tag(NAME_TAG_NAME, NAME_TAG_SIZE,
                    moab::MB_TAG_SPARSE,
                    moab::MB_TYPE_OPAQUE,
                    NULL, false);

  facetingTolTag = get_tag(FACETING_TOL_TAG_NAME, 1,
                           moab::MB_TAG_SPARSE, moab::MB_TYPE_DOUBLE);

  // search for a tag that has the faceting tolerance
  Range tagged_sets;
  double facet_tol_tagvalue = 0;
  bool other_set_tagged = false, root_tagged = false;

  // get list of entity sets that are tagged with faceting tolerance
  // (possibly empty set)
  rval = ErrorCode(MBI->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &facetingTolTag,
                                                     NULL, 1, tagged_sets));
  // if NOT empty set
  if (DAG_SUCCESS == rval && !tagged_sets.empty()) {
    rval = ErrorCode(MBI->tag_get_data(facetingTolTag, &(*tagged_sets.begin()),
                                       1, &facet_tol_tagvalue));
    if (DAG_SUCCESS != rval)
      return rval;
    other_set_tagged = true;
  } else if (DAG_SUCCESS == rval) {
    // check to see if interface is tagged
    EntityHandle root = 0;
    rval = ErrorCode(MBI->tag_get_data(facetingTolTag, &root,
                                       1, &facet_tol_tagvalue));
    if (DAG_SUCCESS == rval)
      root_tagged = true;
    else
      rval = DAG_SUCCESS;
  }

  if ((root_tagged || other_set_tagged) && facet_tol_tagvalue > 0) {
    facetingTolerance = facet_tol_tagvalue;
  }

  // initialize ray_tracer
  std::cout << "Initializing the GeomQueryTool..." << std::endl;
  errHandler->checkSetErr(GTT->find_geomsets(),
                          "Failed to find the geometry sets");

  std::cout << "Using faceting tolerance: " << facetingTolerance << std::endl;

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
    return ErrorCode (ray_tracer->ray_fire(volume, point, dir, next_surf, next_surf_dist,
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

// *****************************************************************************
// SECTION III
// *****************************************************************************

EntityHandle DagMCmoab::entity_by_id(int dimension, int id) {
  return GTT->entity_by_id(dimension, id);
}

int DagMCmoab::id_by_index(int dimension, int index) {
  EntityHandle h = entity_by_index(dimension, index);
  if (!h)
    return 0;

  int result = 0;
  MBI->tag_get_data(GTT->get_gid_tag(), &h, 1, &result);
  return result;
}

int DagMCmoab::get_entity_id(EntityHandle this_ent) {
  return GTT->global_id(this_ent);
}

ErrorCode DagMCmoab::build_indices(Range& surfs, Range& vols) {
  ErrorCode rval = DAG_SUCCESS;


  if (surfs.size() == 0 || vols.size() == 0) {
    std::cout << "Volumes or Surfaces not found" << std::endl;
    return  DAG_ENTITY_NOT_FOUND;
  }
  setOffset = std::min(*surfs.begin(), *vols.begin());
  // surf/vol offsets are just first handles
  EntityHandle tmp_offset = std::max(surfs.back(), vols.back());

  // set size
  entIndices.resize(tmp_offset - setOffset + 1);

  // store surf/vol handles lists (surf/vol by index) and
  // index by handle lists
  surf_handles().resize(surfs.size() + 1);
  std::vector<EntityHandle>::iterator iter = surf_handles().begin();
  // MCNP wants a 1-based index but C++ has a 0-based index. So we need to set
  // the first value to 0 and then start at the next position in the vector
  // (iter++) thereafter.
  *(iter++) = 0;
  std::copy(surfs.begin(), surfs.end(), iter);
  int idx = 1;
  for (Range::iterator rit = surfs.begin(); rit != surfs.end(); ++rit)
    entIndices[*rit - setOffset] = idx++;

  vol_handles().resize(vols.size() + 1);
  iter = vol_handles().begin();

  // MCNP wants a 1-based index but C++ has a 0-based index. So we need to set
  // the first value to 0 and then start at the next position in the vector
  // (iter++) thereafter.
  *(iter++) = 0;
  std::copy(vols.begin(), vols.end(), iter);
  idx = 1;
  for (Range::iterator rit = vols.begin(); rit != vols.end(); ++rit)
    entIndices[*rit - setOffset] = idx++;

  // get group handles
  Tag category_tag = get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE,
                             moab::MB_TAG_SPARSE, moab::MB_TYPE_OPAQUE);
  char group_category[CATEGORY_TAG_SIZE];
  std::fill(group_category, group_category + CATEGORY_TAG_SIZE, '\0');
  sprintf(group_category, "%s", "Group");
  const void* const group_val[] = {&group_category};
  Range groups;
  rval = ErrorCode(MBI->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &category_tag,
                                                     group_val, 1, groups));
  if (DAG_SUCCESS != rval)
    return rval;
  group_handles().resize(groups.size() + 1);
  group_handles()[0] = 0;
  std::copy(groups.begin(), groups.end(), &group_handles()[1]);

  return DAG_SUCCESS;
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

ErrorCode DagMCmoab::write_mesh(const char* ffile,
                                const int flen) {
  ErrorCode rval;

  // write out a mesh file if requested
  if (ffile && 0 < flen) {
    rval = ErrorCode(MBI->write_mesh(ffile));
    if (DAG_SUCCESS != rval) {
      std::cerr << "Failed to write mesh to " << ffile << "." << std::endl;
      return rval;
    }
  }

  return DAG_SUCCESS;
}

// *****************************************************************************
// SECTION V: Metadata handling
// *****************************************************************************

ErrorCode DagMCmoab::get_group_name(EntityHandle group_set, std::string& name) {
  ErrorCode rval;
  const void* v = NULL;
  int ignored;
  rval = ErrorCode(MBI->tag_get_by_ptr(name_tag(), &group_set, 1, &v, &ignored));
  if (DAG_SUCCESS != rval)
    return rval;
  name = static_cast<const char*>(v);
  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::append_packed_string(Tag tag, EntityHandle eh,
                                          std::string& new_string) {
  // When properties have multiple values, the values are tagged in a single character array
  // with the different values separated by null characters
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = ErrorCode(MBI->tag_get_by_ptr(tag, &eh, 1, &p, &len));
  if (rval == DAG_TAG_NOT_FOUND) {
    // This is the first entry, and can be set directly
    p = new_string.c_str();
    return ErrorCode(MBI->tag_clear_data(tag, &eh, 1, p, new_string.length() + 1));
  } else if (rval != DAG_SUCCESS)
    return rval;
  else {
    str = static_cast<const char*>(p);
  }

  // append a new value for the property to the existing property string
  unsigned int tail_len = new_string.length() + null_delimiter_length;
  int new_len = tail_len + len;

  char* new_packed_string = new char[ new_len ];
  memcpy(new_packed_string, str, len);
  memcpy(new_packed_string + len, new_string.c_str(), tail_len);

  p = new_packed_string;
  rval = ErrorCode(MBI->tag_set_by_ptr(tag, &eh, 1, &p, &new_len));
  delete[] new_packed_string;
  return rval;
}

ErrorCode DagMCmoab::unpack_packed_string(Tag tag, EntityHandle eh,
                                          std::vector< std::string >& values) {
  ErrorCode rval;
  const void* p;
  const char* str;
  int len;
  rval = ErrorCode(MBI->tag_get_by_ptr(tag, &eh, 1, &p, &len));
  if (rval != DAG_SUCCESS)
    return rval;
  str = static_cast<const char*>(p);
  int idx = 0;
  while (idx < len) {
    std::string item(str + idx);
    values.push_back(item);
    idx += item.length() + null_delimiter_length;
  }
  return DAG_SUCCESS;
}

ErrorCode DagMCmoab::parse_properties(const std::vector<std::string>& keywords,
                                      const std::map<std::string, std::string>& keyword_synonyms,
                                      const char* delimiters) {
  ErrorCode rval;

  // master keyword map, mapping user-set words in cubit to canonical property names
  std::map< std::string, std::string > keyword_map(keyword_synonyms);

  for (std::vector<std::string>::const_iterator i = keywords.begin();
       i != keywords.end(); ++i) {
    keyword_map[*i] = *i;
  }

  // the set of all canonical property names
  std::set< std::string > prop_names;
  for (prop_map::iterator i = keyword_map.begin();
       i != keyword_map.end(); ++i) {
    prop_names.insert((*i).second);
  }

  // set up DagMC's property tags based on what's been requested
  for (std::set<std::string>::iterator i = prop_names.begin();
       i != prop_names.end(); ++i) {
    std::string tagname("DAGMCPROP_");
    tagname += (*i);

    Tag new_tag;
    rval = ErrorCode(MBI->tag_get_handle(tagname.c_str(), 0,
                                         moab::MB_TYPE_OPAQUE, new_tag,
                                         (moab::MB_TAG_SPARSE |
                                          moab::MB_TAG_VARLEN |
                                          moab::MB_TAG_CREAT)));
    if (DAG_SUCCESS != rval)
      return rval;
    property_tagmap[(*i)] = new_tag;
  }

  // now that the keywords and tags are ready, iterate over all the actual geometry groups
  for (std::vector<EntityHandle>::iterator grp = group_handles().begin();
       grp != group_handles().end(); ++grp) {

    prop_map properties;
    rval = parse_group_name(*grp, properties, delimiters);
    if (rval == DAG_TAG_NOT_FOUND)
      continue;
    else if (rval != DAG_SUCCESS)
      return rval;

    Range grp_sets;
    rval = ErrorCode(MBI->get_entities_by_type(*grp, moab::MBENTITYSET, grp_sets));
    if (DAG_SUCCESS != rval)
      return rval;
    if (grp_sets.size() == 0)
      continue;

    for (prop_map::iterator i = properties.begin();
         i != properties.end(); ++i) {
      std::string groupkey = (*i).first;
      std::string groupval = (*i).second;

      if (property_tagmap.find(groupkey) != property_tagmap.end()) {
        Tag proptag = property_tagmap[groupkey];
        const unsigned int groupsize = grp_sets.size();
        for (unsigned int j = 0; j < groupsize; ++j) {
          rval = append_packed_string(proptag, grp_sets[j], groupval);
          if (DAG_SUCCESS != rval)
            return rval;
        }
      }
    }
  }
  return DAG_SUCCESS;
}

/* ErrorCode DagMC::prop_value(EntityHandle eh, const std::string& prop, std::string& value) {
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return DAG_TAG_NOT_FOUND;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = MBI->tag_get_by_ptr(proptag, &eh, 1, &data, &ignored);
  if (rval != DAG_SUCCESS)
    return rval;
  value = static_cast<const char*>(data);
  return DAG_SUCCESS;
} */

bool DagMCmoab::has_prop(EntityHandle eh, const std::string& prop) {
  ErrorCode rval;

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return false;
  }

  Tag proptag = (*it).second;
  const void* data;
  int ignored;

  rval = ErrorCode(MBI->tag_get_by_ptr(proptag, &eh, 1, &data, &ignored));
  return (rval == DAG_SUCCESS);

}

/* ErrorCode DagMC::get_all_prop_values(const std::string& prop, std::vector<std::string>& return_list) {
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;

  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &proptag, NULL, 1, all_ents);
  if (DAG_SUCCESS != rval)
    return rval;

  std::set<std::string> unique_values;
  for (Range::iterator i = all_ents.begin(); i != all_ents.end(); ++i) {
    std::vector<std::string> values;
    rval = prop_values(*i, prop, values);
    if (DAG_SUCCESS != rval)
      return rval;
    unique_values.insert(values.begin(), values.end());
  }

  return_list.assign(unique_values.begin(), unique_values.end());
  return DAG_SUCCESS;
} */

/* ErrorCode DagMC::entities_by_property(const std::string& prop, std::vector<EntityHandle>& return_list,
                                      int dimension, const std::string* value) {
  ErrorCode rval;
  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return MB_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;
  Range all_ents;

  // Note that we cannot specify values for proptag here-- the passed value,
  // if it exists, may be only a subset of the packed string representation
  // of this tag.
  Tag tags[2] = {proptag, GTT->get_geom_tag()};
  void* vals[2] = {NULL, (dimension != 0) ? &dimension : NULL };
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, tags, vals, 2, all_ents);
  if (DAG_SUCCESS != rval)
    return rval;

  std::set<EntityHandle> handles;
  for (Range::iterator i = all_ents.begin(); i != all_ents.end(); ++i) {
    std::vector<std::string> values;
    rval = prop_values(*i, prop, values);
    if (DAG_SUCCESS != rval)
      return rval;
    if (value) {
      if (std::find(values.begin(), values.end(), *value) != values.end()) {
        handles.insert(*i);
      }
    } else {
      handles.insert(*i);
    }
  }

  return_list.assign(handles.begin(), handles.end());
  return DAG_SUCCESS;
} */

bool DagMCmoab::is_implicit_complement(EntityHandle volume) {
  return GTT->is_implicit_complement(volume);
}

Tag DagMCmoab::get_tag(const char* name, int size, TagType store,
                       DataType type, const void* def_value,
                       bool create_if_missing) {
  Tag retval = 0;
  unsigned flags = store | moab::MB_TAG_CREAT;
  // NOTE: this function seems to be broken in that create_if_missing has
  // the opposite meaning from what its name implies.  However, changing the
  // behavior causes tests to fail, so I'm leaving the existing behavior
  // in place.  -- j.kraftcheck.
  if (!create_if_missing)
    flags |= moab::MB_TAG_EXCL;
  ErrorCode result = ErrorCode(MBI->tag_get_handle(name, size, type, retval, flags, def_value));
  if (create_if_missing && DAG_SUCCESS != result)
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;

  return retval;
}


} // namespace DAGMC
