#ifndef DAGMC_MOAB_HPP
#define DAGMC_MOAB_HPP

// Headers
#include "DagMCBase.hpp"
#include "MoabInterface.hpp"
#include "DefaultRayTracer.hpp"

//Forward  class declarations
class RefEntity;

// What is this for?
struct DagmcVolData {
  int mat_id;
  double density, importance;
  std::string comp_name;
};

namespace DAGMC {

class DagMCmoab : public DagMCBase {

 public:

  // Constructor
  DagMCmoab(std::shared_ptr<moab::Interface> mb_impl = nullptr, double overlap_tolerance = 0., double numerical_precision = .001);

  // Deprecated Constructor
  [[ deprecated("Replaced by DagMC(std::shared_ptr<Interface> mb_impl, ... )") ]]
  DagMCmoab(moab::Interface* mb_impl, double overlap_tolerance = 0., double numerical_precision = .001);

  // Destructor
  ~DagMCmoab() {};

  /** Get subversion revision of this file (DagMCmoab.hpp) */
  [[deprecated]]
  static unsigned int interface_revision() { return 0; }

  /** Git revision of DAGMC */
  inline std::string git_sha() { return DAGMC_GIT_SHA; }

  void init(double overlap_tolerance, double numerical_precision);

  /* SECTION III: Indexing & Cross-referencing */

  /* Most calling apps refer to geometric entities with a combination of
   *  base-1/0 ordinal index (or rank) and global ID (or name).
   *  DagMC also has an internal EntityHandle reference to each geometric entity.
   *  These method provide ways to translate from one to the other.
   */

  /** map from dimension & global ID to EntityHandle */
  EntityHandle entity_by_id(int dimension, int id) override;
  /** map from dimension & base-1 ordinal index to EntityHandle */
  EntityHandle entity_by_index(int dimension, int index) override;
  /** map from dimension & base-1 ordinal index to global ID */
  int id_by_index(int dimension, int index) override;
  /** PPHW: Missing dim & global ID ==> base-1 ordinal index */
  /** map from EntityHandle to base-1 ordinal index */
  int index_by_handle(EntityHandle handle) override;
  /** map from EntityHandle to global ID */
  int get_entity_id(EntityHandle this_ent) override;
  /** \brief get number of geometric sets corresponding to geometry of specified dimension
   *
   * For a given dimension (e.g. dimension=3 for volumes, dimension=2 for surfaces)
   * return the number of entities of that dimension
   *\param dimension the dimensionality of the entities in question
  *\return integer number of entities of that dimension
  */
  unsigned int num_entities(int dimension) override;

  /* SECTION IV: Handling DagMC settings */

  /** retrieve faceting tolerance */
  double faceting_tolerance();

  /* SECTION V: Metadata handling */

  /** \brief Parse properties from group names per metadata syntax standard
   *
   *  @param keywords A list of keywords to parse.  These are considered the canonical
   *                  names of the properties, and constitute the valid inputs to
   *                  has_prop() and prop_value().
   *  @param delimiters An array of characters the routine will use to split the groupname
   *                    into properties.
   *  @param synonyms An optional mapping of synonym keywords to canonical keywords.
   *                  This allows more than one group name keyword to take on the same
   *                  meaning
   *                  e.g. if synonyms["rest.of.world"] = "graveyard", then volumes
   *                  in the "rest.of.world" group will behave as if they were in a
   *                  group named "graveyard".
   */
  ErrorCode parse_properties(const std::vector<std::string>& keywords,
                             const std::map<std::string, std::string>& synonyms = no_synonyms,
                             const char* delimiters = "_") override;

  /** Get the value of a property on a volume or surface
   *
   *  @param eh The entity handle to get a property value on
   *  @param prop The canonical property name
   *  @param value Output parameter, the value of the property.  If no value was
   *               set on the handle, this will be the empty string.
   *  @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *          MOAB, or MB_SUCCESS if successful
   */
  ErrorCode prop_value(EntityHandle eh, const std::string& prop, std::string& value) override;

  /** \brief Get the value of a property on a volume or surface
   *
   *  @param eh The entity handle to get a property value on
   *  @param prop The canonical property name
   *  @param value Output parameter, the value of the property.  If no value was
   *               set on the handle, this will be the empty string.
   *  @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *          MOAB, or MB_SUCCESS if successful
   */
  bool has_prop(EntityHandle eh, const std::string& prop) override;

  /** Detect all the property keywords that appear in the loaded geometry
  *
  *  @param keywords_out The result list of keywords.  This list could be
  *        validly passed to parse_properties().
  */
  ErrorCode detect_available_props(std::vector<std::string>& keywords_out, const char* delimiters = "_")  override;

  /** Get the value of a property on a volume or surface
   *
   *  @param eh The entity handle to get a property value on
   *  @param prop The canonical property name
   *  @param value Output parameter, the value of the property.  If no value was
   *               set on the handle, this will be the empty string.
   *  @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *          MOAB, or MB_SUCCESS if successful
   */
  ErrorCode prop_values(EntityHandle eh, const std::string& prop,
                        std::vector<std::string>& value) override;

  // Non-inherited public!

  // To-Do is this needed? Find where called.
  /** get the tag for the "name" of a surface == global ID */
  moab::Tag name_tag() { return moab_interface->name_tag(); }

  /** Get the tag used to associate OBB trees with geometry in load_file(..).
   * not sure what to do about the obb_tag, GTT has no concept of an obb_tag on EntitySets - PCS
   */
  moab::Tag obb_tag() { return NULL; }
  moab::Tag geom_tag() { return (geom_tool())->get_geom_tag(); }
  moab::Tag id_tag() { return (geom_tool())->get_gid_tag(); }
  moab::Tag sense_tag() { return (geom_tool())->get_sense_tag(); }

 private:

  /** Return true if a volume or surface has the named property set upon it
       *
       *  @param eh The entity handle to query
       *  @param prop The canonical property name
       *  @retrun True if the handle has the property set, or false if not.
       *          False is also returned if a MOAB error occurs.
   */

  ErrorCode get_all_prop_values(const std::string& prop, std::vector<std::string>& return_list);

  /** Get a list of all entities which have a given property
   *
   *  @param prop The canonical property name
   *  @param return_list Output param, a list of entity handles that have this property
   *  @param dimension If nonzero, entities returned will be restricted to the given dimension,
   *                   i.e. 2 for surfaces and 3 for volumes
   *  @parm value If non-NULL, only entities for which the property takes on this value will be returned.
   *  @return MB_TAG_NOT_FOUND if prop is invalid.  Otherwise return any errors from
   *          MOAB, or MB_SUCCESS if successful
   */
  ErrorCode entities_by_property(const std::string& prop, std::vector<EntityHandle>& return_list,
                                 int dimension = 0, const std::string* value = NULL);

  // ***************************************************************************
  // SECTION VI: Other
  // ***************************************************************************

 public:

  /** get the corners of the OBB for a given volume */
  ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]) override;

  /** get the center point and three vectors for the OBB of a given volume */
  ErrorCode getobb(EntityHandle volume, double center[3],
                   double axis1[3], double axis2[3], double axis3[3]) override;

  /** get the root of the obbtree for a given entity */
  ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle& root);

  // Public non-inherited: MOAB  getter methods

  std::shared_ptr<moab::GeomTopoTool> geom_tool() {return moab_interface->gtt();}

  moab::OrientedBoxTreeTool* obb_tree() {return (geom_tool())->obb_tree();}

  /** Get the instance of MOAB used by functions in this file. */
  moab::Interface* moab_instance() { return moab_interface->mesh_ptr();}

  std::shared_ptr<moab::Interface> moab_instance_sptr() {
    return moab_interface->mesh_sptr();
  }

  // ***************************************************************************
  // MEMBER DATA
  // ***************************************************************************

 private:

  // Interface to MOAB
  std::shared_ptr<MoabInterface> moab_interface;

  /** corresponding geometric entities; also indexed like rootSets */
  std::vector<RefEntity*> geomEntities;

  /* metadata */
  char implComplName[NAME_TAG_SIZE];

  /** vectors for point_in_volume: */
  std::vector<double> disList;
  std::vector<int> dirList;
  std::vector<EntityHandle> surList, facList;

}; // End DagMC class definition


// ***************************************************************************
// SECTION VI inline functions
// ***************************************************************************

inline ErrorCode DagMCmoab::getobb(EntityHandle volume, double minPt[3], double maxPt[3]) {
  errHandler->checkSetErr((geom_tool())->get_bounding_coords(volume, minPt, maxPt),
                          "Failed to get obb for volume");
  return DAG_SUCCESS;
}

inline ErrorCode DagMCmoab::getobb(EntityHandle volume, double center[3],
                                   double axis1[3], double axis2[3], double axis3[3]) {
  errHandler->checkSetErr((geom_tool())->get_obb(volume, center, axis1, axis2, axis3),
                          "Failed to get obb for volume");
  return DAG_SUCCESS;
}

inline ErrorCode DagMCmoab::get_root(EntityHandle vol_or_surf, EntityHandle& root) {
  errHandler->checkSetErr((geom_tool())->get_root(vol_or_surf, root),
                          "Failed to get obb root set of volume or surface");
  return DAG_SUCCESS;
}

} // namespace DAGMC

#endif
