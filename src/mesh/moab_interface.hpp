#ifndef DAG_MOAB_INTERFACE_HPP
#define DAG_MOAB_INTERFACE_HPP

#include "mesh_interface.hpp"
#include "Moab.hpp"
#include "Error.hpp"

namespace DAGMC {

class ExternalMOAB : public MeshContainer<moab::Interface> {

 public:
  ExternalMOAB(moab::Interface* moabPtrIn) { moabPtr = moabPtrIn; };
  ~ExternalMOAB() {};
  moab::Interface& mesh() override { return *moabPtr; };
  const moab::Interface& const_mesh() const override { return *moabPtr; };
  moab::Interface* mesh_ptr() { return moabPtr; };
  bool isNull() { return (moabPtr == nullptr); };

 private:
  moab::Interface* moabPtr;

};

class ExternalSharedMOAB : public MeshContainer<moab::Interface> {

 public:
  ExternalSharedMOAB(std::shared_ptr<Interface> MBI_shared_ptr) {
    moabPtr = MBI_shared_ptr;
  }
  ~ExternalSharedMOAB() {};

  moab::Interface& mesh()  override { return *moabPtr; };
  const moab::Interface& const_mesh() const override { return *moabPtr; };
  moab::Interface* ptr() override { return moabPtr.get(); };
  std::shared_ptr<moab::Interface> sptr() override { return moabPtr; };
  bool isNull() override { return (moabPtr == nullptr); };

 private:
  std::shared_ptr<moab::Interface> moabPtr;

};


class InternalMOAB : public MeshContainer<moab::Interface> {

 public:
  InternalMOAB() {
    moabPtr = std::make_shared<moab::Core>();
  };
  ~InternalMOAB() {
    moabPtr->delete_mesh();
  };

  moab::Interface& mesh()  override { return *moabPtr; };
  const moab::Interface& const_mesh() const override { return *moabPtr; };
  moab::Interface* ptr() override { return moabPtr.get(); };
  std::shared_ptr<moab::Interface> sptr() override { return moabPtr; };
  bool isNull() override { return (moabPtr == nullptr); };

 private:
  std::shared_ptr<moab::Interface> moabPtr;

};


class MoabInterface : public MeshInterface {

 public:
  MoabInterface(moab::Interface* moabPtrIn);
  MoabInterface(std::shared_ptr<moab::Interface> moabSharedPtrIn);
  ~MoabInterface() {};

  void init();

  // Read and write to file
  bool load(std::string filename) override;
  bool write(std::string filename) override;

  // Finish the setup of the geometry from an open file
  bool setup_geom();
  bool setup_indices();

  // Settings
  bool set_faceting_tol();
  double get_faceting_tol() { return facetingTolerance; };

  // Methods for retrieving and updating metadata
  bool update_properties(std::set< std::string >& prop_names, const char* delimiters);
  bool get_keywords(std::vector<std::string>& keywords_list, const char* delimiters);
  bool get_property(EntityHandle eh, const std::string& prop, std::string& value);
  bool get_properties(EntityHandle eh, const std::string& prop, std::vector< std::string >& values);
  bool has_property(EntityHandle eh, const std::string& prop);
  bool get_ents_and_vals_with_prop(const std::string& prop,
                                   std::set<EntityHandle>& handles,
                                   std::set<std::string>& unique_vals,
                                   bool checkval = false,
                                   int dimension = 0,
                                   std::string value = "");

  // Retrieve references to moab
  moab::Interface& moab() { return container->mesh(); };
  moab::Interface* moab_ptr() { return container->ptr(); };
  std::shared_ptr<moab::Interface> moab_sptr() { return container->sptr(); };

  // Return a copy of the topo tool
  std::shared_ptr<GeomTopoTool> gtt() { return GTT; };

  // Return error code, but cast as DAGMC error code.
  ErrorCode code() { return ErrorCode(rval); };

  // Get the name tag
  // To-do: what is this for?
  Tag name_tag() { return nameTag; };

  // Indexing methods

  /** map from dimension & global ID to EntityHandle */
  EntityHandle entity_by_id(int dimension, int id);
  /** map from dimension & base-1 ordinal index to EntityHandle */
  EntityHandle entity_by_index(int dimension, int index);

  /** map from dimension & base-1 ordinal index to global ID */
  int id_by_index(int dimension, int index);
  /** PPHW: Missing dim & global ID ==> base-1 ordinal index */
  /** map from EntityHandle to base-1 ordinal index */
  int index_by_handle(EntityHandle handle);
  /** map from EntityHandle to global ID */
  int get_entity_id(EntityHandle this_ent);


  /** \brief get number of geometric sets corresponding to geometry of specified dimension
   *
   * For a given dimension (e.g. dimension=3 for volumes, dimension=2 for surfaces)
   * return the number of entities of that dimension
   *\param dimension the dimensionality of the entities in question
  *\return integer number of entities of that dimension
  */
  unsigned int num_entities(int dimension);


 private:

  /** a common type within the property and group name functions */
  typedef std::map<std::string, std::string> prop_map;

  bool build_indices(Range& surfs, Range& vols);

  // Return a map of tags given property names
  bool set_tagmap(std::set< std::string >& prop_names);

  /** \brief Parse a group name into a set of key:value pairs */
  bool parse_group_name(EntityHandle group_set, prop_map& result, const char* delimiters = "_");

  /** \brief tokenize the metadata stored in group names
  * - basically borrowed from ReadCGM.cpp.
  * Called by parse_group_name
  * */
  void tokenize(const std::string& str, std::vector<std::string>& tokens,
                const char* delimiters = "_") const;

  /** Add a string value to a property tag for a given entity */
  bool append_packed_string(Tag, EntityHandle, std::string&);

  // Methods for setting metadata
  bool set_tag(Tag tag, EntityHandle eh, std::string& new_string);
  bool set_tag_data(Tag tag, const EntityHandle* entityPtr,
                    int num_handles, const void* const tag_data, int len);
  bool append_group_properties(const char* delimiters);


  // Methods for fetching metadata values
  bool get_group_handles(std::vector<EntityHandle>& group_handles);
  bool get_group_props(EntityHandle group, std::string& name);
  bool get_entity_sets(EntityHandle group, Range& group_sets);
  bool get_tag(std::string& tagname, Tag& tag);
  bool get_tag_data(const Tag& tag, const EntityHandle* entityPtr,
                    const int num_handles, void* tag_data);
  bool get_tag_data_vec(Tag tag, EntityHandle eh, std::vector<std::string>& values);
  bool get_tag_data_arr(const Tag& tag, const EntityHandle* entityPtr,
                        const int num_handles, const void** tag_data, int* len);
  bool get_tag_name(Tag tag, EntityHandle eh, std::string& name);
  Tag  get_tag(const char* name, int size, TagType store, DataType type,
               const void* def_value = NULL, bool create_if_missing = true);
  bool get_prop_tag(const std::string& prop, Tag& proptag);
  bool get_tagged_entity_sets(EntityHandle group, Tag tag, Range& group_sets);
  bool get_tagged_entity_sets(EntityHandle group, std::vector<Tag> tags,
                              const void* const* vals, Range& group_sets);

  // Convenience methods for accessing entity handles
  std::vector<EntityHandle>& surf_handles() {
    return entHandles[surfs_handle_idx];
  };
  std::vector<EntityHandle>& vol_handles() {
    return entHandles[vols_handle_idx];
  };
  std::vector<EntityHandle>& group_handles() {
    return entHandles[groups_handle_idx];
  };

  // Reset the state of the MOAB error code.
  void reset_code() { rval = moab::MB_SUCCESS; };

  // Container for the mesh
  std::shared_ptr<MeshContainer<moab::Interface> > container;

  // Pointer to an instance of a GeomTopoTool
  std::shared_ptr<GeomTopoTool> GTT;

  // Store MOAB return values
  moab::ErrorCode rval;

  Tag nameTag;
  Tag facetingTolTag;

  static const int null_delimiter_length = 1;

  static const int vertex_handle_idx = 0;
  static const int curve_handle_idx = 1;
  static const int surfs_handle_idx = 2;
  static const int vols_handle_idx = 3;
  static const int groups_handle_idx = 4;

  double facetingTolerance;

  /** store some lists indexed by handle */
  std::vector<EntityHandle> entHandles[5];

  /** lowest-valued handle among entity sets representing surfs and vols */
  EntityHandle setOffset;

  /** entity index (contiguous 1-N indices); indexed like rootSets */
  std::vector<int> entIndices;

  /** map from the canonical property names to the tags representing them */
  std::map<std::string, Tag> property_tagmap;

};


inline EntityHandle MoabInterface::entity_by_index(int dimension, int index) {
  assert(2 <= dimension && 3 >= dimension && (unsigned)index < entHandles[dimension].size());
  return entHandles[dimension][index];
}

inline int MoabInterface::index_by_handle(EntityHandle handle) {
  assert(handle - setOffset < entIndices.size());
  return entIndices[handle - setOffset];
}

inline unsigned int MoabInterface::num_entities(int dimension) {
  assert(vertex_handle_idx <= dimension && groups_handle_idx >= dimension);
  return entHandles[dimension].size() - 1;
}

}


#endif
