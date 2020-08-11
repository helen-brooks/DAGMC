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

// initialise the obb tree
ErrorCode DagMCBase::init_OBBTree() {

  // find all geometry sets
  errHandler->checkSetErr(find_geomsets(),
                          "Could not find the geometry sets");

  // implicit complement
  errHandler->checkSetErr(setup_impl_compl(),
                          "Failed to setup the implicit compliment");

  // build obbs
  errHandler->checkSetErr(setup_obbs(), "Failed to setup the OBBs");

  // setup indices
  errHandler->checkSetErr(setup_indices(), "Failed to setup problem indices");

  return DAG_SUCCESS;
}

// helper function to load the existing contents of a MOAB instance into DAGMC
ErrorCode DagMCBase::load_existing_contents() {
  return finish_loading();
}

// setups of the indices for the problem, builds a list of surface and volumes
// indices
ErrorCode DagMCBase::setup_indices() {
  Range surfs, vols;
  errHandler->checkSetErr(setup_geometry(surfs, vols),
                          "Failed to setup geometry");
  // build the various index vectors used for efficiency
  errHandler->checkSetErr(build_indices(surfs, vols),
                          "Failed to build surface/volume indices");

  return errHandler->code();
}

ErrorCode DagMCBase::detect_available_props(std::vector<std::string>& keywords_list,
                                            const char* delimiters) {
  ErrorCode rval;
  std::set< std::string > keywords;
  for (std::vector<EntityHandle>::const_iterator grp = group_handles().begin();
       grp != group_handles().end(); ++grp) {
    std::map< std::string, std::string > properties;
    rval = parse_group_name(*grp, properties, delimiters);
    if (rval == DAG_TAG_NOT_FOUND)
      continue;
    else if (rval != DAG_SUCCESS)
      return rval;

    for (prop_map::iterator i = properties.begin();
         i != properties.end(); ++i) {
      keywords.insert((*i).first);
    }
  }
  keywords_list.assign(keywords.begin(), keywords.end());
  return DAG_SUCCESS;
}

ErrorCode DagMCBase::prop_values(EntityHandle eh, const std::string& prop,
                                 std::vector< std::string >& values) {

  std::map<std::string, Tag>::iterator it = property_tagmap.find(prop);
  if (it == property_tagmap.end()) {
    return DAG_TAG_NOT_FOUND;
  }
  Tag proptag = (*it).second;

  return unpack_packed_string(proptag, eh, values);

}

// *****************************************************************************
// Protected methods
// *****************************************************************************

ErrorCode DagMCBase::parse_group_name(EntityHandle group_set, prop_map& result,
                                      const char* delimiters) {
  ErrorCode rval;
  std::string group_name;
  rval = get_group_name(group_set, group_name);
  if (rval != DAG_SUCCESS)
    return rval;

  std::vector< std::string > group_tokens;
  tokenize(group_name, group_tokens, delimiters);

  // iterate over all the keyword positions
  // keywords are even indices, their values (optional) are odd indices
  for (unsigned int i = 0; i < group_tokens.size(); i += 2) {
    std::string groupkey = group_tokens[i];
    std::string groupval;
    if (i < group_tokens.size() - 1)
      groupval = group_tokens[i + 1];
    result[groupkey] = groupval;
  }
  return DAG_SUCCESS;
}

void DagMCBase::tokenize(const std::string& str,
                         std::vector<std::string>& tokens,
                         const char* delimiters) const {
  std::string::size_type last = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos  = str.find_first_of(delimiters, last);
  if (std::string::npos == pos)
    tokens.push_back(str);
  else
    while (std::string::npos != pos && std::string::npos != last) {
      tokens.push_back(str.substr(last, pos - last));
      last = str.find_first_not_of(delimiters, pos);
      pos  = str.find_first_of(delimiters, last);
      if (std::string::npos == pos)
        pos = str.size();
    }
}