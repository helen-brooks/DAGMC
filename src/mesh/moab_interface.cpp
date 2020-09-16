#include "moab_interface.hpp"

namespace DAGMC {

#define FACETING_TOL_TAG_NAME "FACETING_TOL"

// *****************************************************************************
// PUBLIC METHODS
// *****************************************************************************

MoabInterface::MoabInterface(moab::Interface* moabPtrIn) :
  null_delimiter_length(1) {
  rval = moab::MB_SUCCESS;
  container = std::make_shared< ExternalMOAB >(moabPtrIn);
}

MoabInterface::MoabInterface(std::shared_ptr<moab::Interface> moabSharedPtrIn) :
  null_delimiter_length(1) {
  rval = moab::MB_SUCCESS;
  if (moabSharedPtrIn != nullptr) {
    container = std::make_shared< ExternalSharedMOAB >(moabSharedPtrIn);
  } else {
    container = std::make_shared< InternalMOAB >();
  }
}

bool MoabInterface::load(std::string filename) {

  std::cout << "Loading file " << filename << std::endl;

  EntityHandle file_set;
  rval = moab().create_meshset(moab::MESHSET_SET, file_set);
  if (moab::MB_SUCCESS != rval)
    return false;

  // load options

  // set entire options string to literal null character
  // Can't we just pass '\0'
  char options[120] = {0};
  rval = moab().load_file(filename.c_str(), &file_set, options, NULL, 0, 0);

  if (moab::MB_UNHANDLED_OPTION == rval) {
    // Some options were unhandled; this is common for loading h5m files.
    // Print a warning if an option was unhandled for a file that does not end in '.h5m'
    std::string ext = ".h5m";
    size_t pos = filename.rfind(ext);
    if (pos == std::string::npos || pos != ext.size()) {
      std::cerr << "DagMC warning: unhandled file loading options." << std::endl;
    }
  } else if (moab::MB_SUCCESS != rval) {
    std::cerr << "DagMC Couldn't read file " << filename << std::endl;

    std::string message;
    if (moab::MB_SUCCESS == moab().get_last_error(message) && !message.empty())
      std::cerr << "Error message: " << message << std::endl;

    return false;
  }

  return true;
}

bool MoabInterface::write(std::string filename) {

  rval = moab().write_mesh(filename.c_str());
  return (rval == moab::MB_SUCCESS);

}

bool MoabInterface::get_faceting_tol(double& facetingTolerance) {

  nameTag = get_tag(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TAG_SPARSE,
                    moab::MB_TYPE_OPAQUE, NULL, false);
  facetingTolTag = get_tag(FACETING_TOL_TAG_NAME, 1, moab::MB_TAG_SPARSE,
                           moab::MB_TYPE_DOUBLE);

  // search for a tag that has the faceting tolerance
  Range tagged_sets;
  double facet_tol_tagvalue = 0;
  bool tagged = false;
  EntityHandle root = 0;

  // get list of entity sets that are tagged with faceting tolerance
  // (possibly empty set)
  if (get_tagged_entity_sets(root, facetingTolTag, tagged_sets)) {

    // Found tagged sets
    if (!tagged_sets.empty()) {
      tagged = get_tag_data(facetingTolTag, &tagged_sets.front(), 1, &facet_tol_tagvalue);
      // Failed to get data
      if (!tagged)
        return false;
    } else {
      // Check to see if interface is tagged
      tagged = get_tag_data(facetingTolTag, &root, 1, &facet_tol_tagvalue);
    }
  }

  // Interface was not tagged
  if (!tagged) {
    // This is fine, just use default value, but need to reset rval
    rval = moab::MB_SUCCESS;
  }
  // Inteface was tagged, but check for positive value
  else if (facet_tol_tagvalue > 0) {
    facetingTolerance = facet_tol_tagvalue;
  }

  std::cout << "Using faceting tolerance: " << facetingTolerance << std::endl;
  return true;
}

bool MoabInterface::get_group_handles(std::vector<EntityHandle>& group_handles) {

  Tag category_tag = get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE,
                             moab::MB_TAG_SPARSE, moab::MB_TYPE_OPAQUE);
  char group_category[CATEGORY_TAG_SIZE];
  std::fill(group_category, group_category + CATEGORY_TAG_SIZE, '\0');
  sprintf(group_category, "%s", "Group");
  const void* const group_val[] = {&group_category};

  Range groups;
  rval = moab().get_entities_by_type_and_tag(0, moab::MBENTITYSET, &category_tag,
                                             group_val, 1, groups);
  if (moab::MB_SUCCESS != rval)
    return false;

  group_handles.resize(groups.size() + 1);
  group_handles[0] = 0;
  std::copy(groups.begin(), groups.end(), &group_handles[1]);

  return true;
}

bool MoabInterface::get_tag(std::string& tagname, Tag& tag) {
  rval = moab().tag_get_handle(tagname.c_str(), 0,
                               moab::MB_TYPE_OPAQUE, tag,
                               (moab::MB_TAG_SPARSE |
                                moab::MB_TAG_VARLEN |
                                moab::MB_TAG_CREAT));
  if (moab::MB_SUCCESS != rval)
    return false;
  else
    return true;

}

bool MoabInterface::get_tag_data(const Tag& tag, const EntityHandle* entityPtr,
                                 const int num_handles, void* tag_data) {

  rval = moab().tag_get_data(tag, entityPtr, num_handles, tag_data);
  if (moab::MB_SUCCESS != rval)
    return false;
  else
    return true;

}

bool MoabInterface::get_tag_data(Tag tag, EntityHandle eh, std::vector<std::string>& values) {

  // //  std::string str;
  // const void* data;
  // int len;
  // rval = moab().tag_get_by_ptr(tag, &eh, 1, &data, &len);

  // if (moab::MB_SUCCESS != rval) return false;

  // const char* c_str = static_cast<const char*>(data);
  // int idx = 0;
  // while (idx < len) {
  //   std::string item(c_str + idx);
  //   values.push_back(item);
  //   idx += item.length() + null_delimiter_length;
  // }
  // return true;


  const void* p;
  const char* str;
  int len;
  rval = moab().tag_get_by_ptr(tag, &eh, 1, &p, &len);
  if (rval != moab::MB_SUCCESS)
    return false;

  str = static_cast<const char*>(p);
  int idx = 0;
  while (idx < len) {
    std::string item(str + idx);
    values.push_back(item);
    idx += item.length() + null_delimiter_length;
  }
  return true;

}

bool MoabInterface::get_tag_name(Tag tag, EntityHandle eh, std::string& name) {

  const void* data;
  int len;
  rval = moab().tag_get_by_ptr(tag, &eh, 1, &data, &len);

  if (moab::MB_SUCCESS != rval)
    return false;

  name = std::string(static_cast<const char*>(data));

  return true;
}

bool MoabInterface::get_group_name(EntityHandle group, std::string& name) {

  return get_tag_name(nameTag, group, name);

}

bool MoabInterface::get_entity_sets(EntityHandle group, Range& group_sets) {

  rval = moab().get_entities_by_type(group, moab::MBENTITYSET, group_sets);
  return (rval == moab::MB_SUCCESS);

}

bool MoabInterface::get_tagged_entity_sets(EntityHandle group, Tag tag, Range& group_sets) {

  return get_tagged_entity_sets(group, std::vector<Tag>(1, tag), group_sets);

}

bool MoabInterface::get_tagged_entity_sets(EntityHandle group, std::vector<Tag> tags, Range& group_sets) {

  rval = moab().get_entities_by_type_and_tag(group, moab::MBENTITYSET, tags.data(), NULL, tags.size(), group_sets);
  return (rval == moab::MB_SUCCESS);
}


bool MoabInterface::set_tag(Tag tag, EntityHandle eh, std::string& new_string) {

  int len = new_string.length() + null_delimiter_length;
  const void* data = new_string.c_str();
  rval = moab().tag_set_by_ptr(tag, &eh, 1, &data, &len);
  return (rval == moab::MB_SUCCESS);

}


// *****************************************************************************
// PRIVATE METHODS
// *****************************************************************************

Tag MoabInterface::get_tag(const char* name, int size, TagType store,
                           DataType type, const void* def_value,
                           bool create_if_missing) {

  unsigned flags = store | moab::MB_TAG_CREAT;
  // NOTE: this function seems to be broken in that create_if_missing has
  // the opposite meaning from what its name implies.  However, changing the
  // behavior causes tests to fail, so I'm leaving the existing behavior
  // in place.  -- j.kraftcheck.
  if (!create_if_missing)
    flags |= moab::MB_TAG_EXCL;

  Tag retval = 0;
  rval = moab().tag_get_handle(name, size, type, retval, flags, def_value);
  if (create_if_missing && moab::MB_SUCCESS != rval)
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;

  return retval;
}


}
