#ifndef DAGMC_MOAB_HPP
#define DAGMC_MOAB_HPP

// Headers
#include "DagMCBase.hpp"
#include "moab_interface.hpp"

//Forward  class declarations
class RefEntity;
class RayTracingInterface;

// What is this for?
struct DagmcVolData {
  int mat_id;
  double density, importance;
  std::string comp_name;
};

namespace DAGMC {

class MoabErrHandler : public ErrorHandler {
 public:
  MoabErrHandler() {};
  ~MoabErrHandler() {};

  void checkSetErr(ErrorCode rval, std::string msg) override {
    _code = ErrorCode(setMoabCode(rval, msg));
    //possible exception handling here?
    return;
  };
  //Wrapper for MOAB macro which contains a return moab::ErrorCode statement
  moab::ErrorCode setMoabCode(ErrorCode rval, std::string msg) {
    moab::ErrorCode mbcode = moab::ErrorCode(rval);
    MB_CHK_SET_ERR(mbcode, msg);
    return mbcode;
  };
};

class DagMCmoab : public DagMCBase {

 public:

  // Constructor
  DagMCmoab(std::shared_ptr<Interface> mb_impl = nullptr, double overlap_tolerance = 0., double numerical_precision = .001);

  // Deprecated Constructor
  [[ deprecated("Replaced by DagMC(std::shared_ptr<Interface> mb_impl, ... )") ]]
  DagMCmoab(Interface* mb_impl, double overlap_tolerance = 0., double numerical_precision = .001);

  // Destructor
  ~DagMCmoab() {};

  /** Get subversion revision of this file (DagMCmoab.hpp) */
  [[deprecated]]
  static unsigned int interface_revision() { return 0; }

  /** Git revision of DAGMC */
  inline std::string git_sha() { return DAGMC_GIT_SHA; }

  /* SECTION I: Geometry Initialization */

  /** \brief Load a geometry description regardless of format
   *
   * This method will load the geometry file with name cfile.
   * In case this is a solid model geometry file, it will pass
   * the facet_tolerance option as guidance for the faceting engine.
   * \param cfile the file name to be loaded
   * \param facet_tolerance the faceting tolerance guidance for the faceting engine
   * \return - MB_SUCCESS if file loads correctly
   *        - other MB ErrorCodes returned from MOAB
   *
   * Note: When loading a prexisting file with an OBB_TREE tag, a number of unspoken
   * things happen that one should be aware of.
   *
   * 1) The file is loaded and when we query the meshset, we find entities with the OBB_TREE tag
   * 2) The OBBTreeTool assumes that any children of the entity being queried in a ray intersect sets
   *      operation are fair game, the surface meshsets have triangles as members, but OBBs as children
   *      but no querying is done, just assumptions that the tags exist.
   */
  ErrorCode load_file(const char* cfile) override;

  /** \brief Use pre-loaded geometry set
   *
   * Works like load_file, but using data that has been externally
   * loaded into DagMC's MOAB instance.
   * Only one of the two functions should be called.
   *
   * TODO: this function should accept a parameter, being the
   * entity set to use for DagMC's data.  Currently DagMC always
   * assumes that all the contents of its MOAB instance belong to it.
   */
  ErrorCode load_existing_contents() override;

  /** \brief finds or creates the implicit complement
   *
   * This method calls the GeomTopoTool->get_implicit_complement which will
   * return the IC if it already exists. If the IC doesn't exist, it will
   * create one.
   */
  ErrorCode setup_impl_compl() override;

  /** \brief initializes the geometry and OBB tree structure for ray firing acceleration
  *
  * This method can be called after load_file to fully initialize DAGMC.  It calls
  * methods to set up the geometry, create the implicit complement, generate an
  * OBB tree from the faceted representation of the geometry, and build the
  * cross-referencing indices.
  */
  ErrorCode init_OBBTree() override;

  /** \brief constructs obb trees for all surfaces and volumes
   *
   * Very thin wrapper around GTT->construct_obb_trees().
   * Constructs obb trees for all surfaces and volumes in the geometry.
   * Called by init_OBBTree
   */
  ErrorCode setup_obbs() override;

  /** Thin wrapper around build_indices().*/
  ErrorCode setup_indices() override;

 private:

  /* SECTION II: Fundamental Geometry Operations/Queries */
  /* The methods in this section are thin wrappers around methods in the
   *  GeometryQueryTool.
   */

  typedef GeomQueryTool::RayHistory RayHistory;

 public:

  ErrorCode ray_fire(const EntityHandle volume, const double ray_start[3],
                     const double ray_dir[3], EntityHandle& next_surf,
                     double& next_surf_dist,
                     RayHistory* history = NULL,
                     double dist_limit = 0, int ray_orientation = 1,
                     OrientedBoxTreeTool::TrvStats* stats = NULL) override;

  ErrorCode point_in_volume(const EntityHandle volume, const double xyz[3],
                            int& result, const double* uvw = NULL,
                            const RayHistory* history = NULL) override;

  ErrorCode point_in_volume_slow(const EntityHandle volume, const double xyz[3],
                                 int& result) override;

  ErrorCode test_volume_boundary(const EntityHandle volume,
                                 const EntityHandle surface,
                                 const double xyz[3], const double uvw[3],
                                 int& result,
                                 const RayHistory* history = NULL) override;

  ErrorCode closest_to_location(EntityHandle volume, const double point[3],
                                double& result, EntityHandle* surface = 0) override;
  ErrorCode measure_volume(EntityHandle volume, double& result) override;

  ErrorCode measure_area(EntityHandle surface, double& result) override;

  ErrorCode surface_sense(EntityHandle volume, int num_surfaces,
                          const EntityHandle* surfaces, int* senses_out) override;

  ErrorCode surface_sense(EntityHandle volume, EntityHandle surface,
                          int& sense_out) override;

  ErrorCode get_angle(EntityHandle surf, const double xyz[3], double angle[3],
                      const RayHistory* history = NULL) override;

  ErrorCode next_vol(EntityHandle surface, EntityHandle old_volume,
                     EntityHandle& new_volume) override;

  /* SECTION III: Indexing & Cross-referencing */
 public:
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
 public:

  /** retrieve faceting tolerance */
  double faceting_tolerance() override {
    return mesh_interface->get_faceting_tol();
  };
  /** retrieve overlap thickness */
  double overlap_thickness() override;
  /** retrieve numerical precision */
  double numerical_precision();
  /** Attempt to set a new overlap thickness tolerance, first checking for sanity */
  void set_overlap_thickness(double new_overlap_thickness) override;

  /** Attempt to set a new numerical precision , first checking for sanity
   *  Use of this function is discouraged; see top of DagMC.cpp
   */
  void set_numerical_precision(double new_precision);


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

  bool is_implicit_complement(EntityHandle volume) override;

  // Non-inherited public!

  // To-Do is this needed? Find where called.
  /** get the tag for the "name" of a surface == global ID */
  Tag name_tag() { return mesh_interface->name_tag(); }

  /** Get the tag used to associate OBB trees with geometry in load_file(..).
   * not sure what to do about the obb_tag, GTT has no concept of an obb_tag on EntitySets - PCS
   */
  Tag obb_tag() { return NULL; }
  Tag geom_tag() { return GTT->get_geom_tag(); }
  Tag id_tag() { return GTT->get_gid_tag(); }
  Tag sense_tag() { return GTT->get_sense_tag(); }

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

  ErrorCode write_mesh(const char* ffile) override;
  /** get the corners of the OBB for a given volume */
  ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]) override;
  /** get the center point and three vectors for the OBB of a given volume */
  ErrorCode getobb(EntityHandle volume, double center[3],
                   double axis1[3], double axis2[3], double axis3[3]) override;
  /** get the root of the obbtree for a given entity */
  ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle& root);

  // public non-inherited: moab getter methods
  OrientedBoxTreeTool* obb_tree() {return GTT->obb_tree();}
  std::shared_ptr<GeomTopoTool> geom_tool() {return GTT;}
  /** Get the instance of MOAB used by functions in this file. */
  Interface* moab_instance() { return mesh_interface->moab_ptr();}
  std::shared_ptr<Interface> moab_instance_sptr() {
    return mesh_interface->moab_sptr();
  }

  // ***************************************************************************
  // MEMBER DATA
  // ***************************************************************************

 private:

  // Interface to MOAB
  std::shared_ptr<MoabInterface> mesh_interface;

  std::shared_ptr<GeomTopoTool> GTT;
  // type alias for ray tracing engine
#ifdef DOUBLE_DOWN
  using RayTracer = RayTracingInterface;
#else
  using RayTracer = GeomQueryTool;
#endif

  std::unique_ptr<RayTracer> ray_tracer;

  /** corresponding geometric entities; also indexed like rootSets */
  std::vector<RefEntity*> geomEntities;

  /* metadata */
  char implComplName[NAME_TAG_SIZE];

  double facetingTolerance;

  /** vectors for point_in_volume: */
  std::vector<double> disList;
  std::vector<int> dirList;
  std::vector<EntityHandle> surList, facList;

}; // End DagMC class definition


// ***************************************************************************
// SECTION VI inline functiosn
// ***************************************************************************

inline ErrorCode DagMCmoab::getobb(EntityHandle volume, double minPt[3], double maxPt[3]) {
  errHandler->checkSetErr(GTT->get_bounding_coords(volume, minPt, maxPt),
                          "Failed to get obb for volume");
  return DAG_SUCCESS;
}

inline ErrorCode DagMCmoab::getobb(EntityHandle volume, double center[3],
                                   double axis1[3], double axis2[3], double axis3[3]) {
  errHandler->checkSetErr(GTT->get_obb(volume, center, axis1, axis2, axis3),
                          "Failed to get obb for volume");
  return DAG_SUCCESS;
}

inline ErrorCode DagMCmoab::get_root(EntityHandle vol_or_surf, EntityHandle& root) {
  errHandler->checkSetErr(GTT->get_root(vol_or_surf, root),
                          "Failed to get obb root set of volume or surface");
  return DAG_SUCCESS;
}

} // namespace DAGMC

#endif
