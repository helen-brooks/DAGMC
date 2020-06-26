//This file was created by H Brooks on 26/06/20
#ifndef DAGMCBASE_HPP

// Header files
#include "moab/Core.hpp"

// Some typedefs for commonly used moab types
// @TODO Generalise and reduce MOAB dependence
typedef moab::ErrorCode ErrorCode;
typedef moab::EntityHandle EntityHandle;

// This is a base class for the DAGMC object whose derivatives will be MOAB- or
// LibMesh- dependent respectively.
class DagMCBase {
 public:
  // Constructor
  DagMCBase();
  // Destructor
  ~DagMCBase();

  // Public methods
  //@TODO Add some pure virtual methods here...

  // SECTION I: Geometry Initialization
  virtual ErrorCode load_file(const char* cfile) = 0;
  virtual ErrorCode setup_impl_compl() = 0;
  virtual ErrorCode setup_indices() = 0;

  // SECTION II: Fundamental Geometry Operations/Queries
  virtual bool is_implicit_complement(EntityHandle volume) = 0;

  // SECTION III: Indexing & Cross-referencing
  virtual EntityHandle entity_by_id(int dimension, int id) = 0;
  virtual EntityHandle entity_by_index(int dimension, int index) = 0;
  virtual int id_by_index(int dimension, int index) = 0;
  /*
    int index_by_handle(EntityHandle handle);
    int get_entity_id(EntityHandle this_ent);
   */
  virtual unsigned int num_entities(int dimension) = 0;

  // SECTION IV: Handling DagMC settings

  // SECTION V: Metadata handling
  virtual ErrorCode parse_properties(const std::vector<std::string>& keywords,
                                     const std::map<std::string, std::string>& synonyms = no_synonyms,
                                     const char* delimiters = "_") = 0;
  virtual ErrorCode prop_values(EntityHandle eh, const std::string& prop,
                                std::vector< std::string >& value) = 0;
  virtual bool has_prop(EntityHandle eh, const std::string& prop) = 0;
}

} // end DAGMC namespace

#define DAGMCBASE_HPP
#endif