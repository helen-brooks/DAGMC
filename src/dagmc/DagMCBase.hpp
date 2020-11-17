//This file was created by H Brooks on 26/06/20
#ifndef DAGMCBASE_HPP
#define DAGMCBASE_HPP

// DAGMC headers
#include "DagMCVersion.hpp"
#include "Common.hpp"
#include "Types.hpp"
#include "MeshInterface.hpp"
#include "RayTracer.hpp"

namespace DAGMC {

/** \brief This is a base class for the DagMC object whose derivatives will be MOAB- or
 * LibMesh- dependent respectively.
 *
 * In section 1, the public interface you will find all the functions needed
 * for problem setup. For the typical MC code, the order of function calls
 * required to fully populate DAGMC ready to run are
 *
 *    1) DAG->load_file();
 *    2) DAG->init_OBBTree();
 *
 * Modifications were made to init_OBBTree which allows the functions of
 * init_OBBTree to be called without having used init_OBBTree. For example
 * if you would like access to be able to call DAG->point_in_volume() but without
 * having an implicit compliment you need only call
 *
 *   1) DAG->load_file();
 *   2) DAG->setup_obbs();
 *
 * Similarly, if you need access to problem indices only, one may call
 * load_file followed by setup_indices.
 *
 *   1) DAG->load_file();
 *   2) DAG->setup_indices();
 */
class DagMCBase {
 public:
  /** Constructor */
  DagMCBase() {};
  /** Destructor */
  ~DagMCBase() {};

  // ***************************************************************************
  // Public methods
  // ***************************************************************************

  /** Return the version of this library */
  static float version(std::string* version_string = NULL);

  // ***************************************************************************
  // SECTION I: Geometry Initialization
  // ***************************************************************************

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
  ErrorCode load_file(const char* cfile);

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
  ErrorCode load_existing_contents();

  /** \brief finds or creates the implicit complement
   *
   * This method calls the GeomTopoTool->get_implicit_complement which will
   * return the IC if it already exists. If the IC doesn't exist, it will
   * create one.
   */
  ErrorCode setup_impl_compl();

  /** \brief initializes the geometry and OBB tree structure for ray firing acceleration
  *
  * This method can be called after load_file to fully initialize DAGMC.  It calls
  * methods to set up the geometry, create the implicit complement, generate an
  * OBB tree from the faceted representation of the geometry, and build the
  * cross-referencing indices.
  */
  ErrorCode init_OBBTree();

  /** \brief constructs obb trees for all surfaces and volumes
   *
   * Thin wrapper around RayTracer::init_obb().
   */
  ErrorCode setup_obbs();

  /** Thin wrapper around MeshInterface::setup_indices().*/
  ErrorCode setup_indices();

  // ***************************************************************************
  // SECTION II: Fundamental Geometry Operations/Queries
  // ***************************************************************************

  ErrorCode ray_fire(const EntityHandle volume, const double ray_start[3],
                     const double ray_dir[3], EntityHandle& next_surf,
                     double& next_surf_dist,
                     moab::GeomQueryTool::RayHistory* history = NULL,
                     double dist_limit = 0, int ray_orientation = 1,
                     moab::OrientedBoxTreeTool::TrvStats* stats = NULL);

  ErrorCode point_in_volume(const EntityHandle volume, const double xyz[3],
                            int& result, const double* uvw = NULL,
                            const moab::GeomQueryTool::RayHistory* history = NULL);

  ErrorCode point_in_volume_slow(const EntityHandle volume, const double xyz[3],
                                 int& result);

  ErrorCode test_volume_boundary(const EntityHandle volume,
                                 const EntityHandle surface,
                                 const double xyz[3], const double uvw[3],
                                 int& result,
                                 const moab::GeomQueryTool::RayHistory* history = NULL);

  ErrorCode closest_to_location(EntityHandle volume, const double point[3],
                                double& result, EntityHandle* surface = 0);
  ErrorCode measure_volume(EntityHandle volume, double& result);

  ErrorCode measure_area(EntityHandle surface, double& result);

  ErrorCode surface_sense(EntityHandle volume, int num_surfaces,
                          const EntityHandle* surfaces, int* senses_out);

  ErrorCode surface_sense(EntityHandle volume, EntityHandle surface,
                          int& sense_out);

  ErrorCode get_angle(EntityHandle surf, const double xyz[3], double angle[3],
                      const moab::GeomQueryTool::RayHistory* history = NULL);

  ErrorCode next_vol(EntityHandle surface, EntityHandle old_volume,
                     EntityHandle& new_volume);

  bool is_implicit_complement(EntityHandle volume);

  // ***************************************************************************
  // SECTION III: Indexing & Cross-referencing
  // ***************************************************************************

  virtual EntityHandle entity_by_id(int dimension, int id) = 0;
  virtual EntityHandle entity_by_index(int dimension, int index) = 0;
  virtual int id_by_index(int dimension, int index) = 0;
  virtual int index_by_handle(EntityHandle handle) = 0;
  virtual int get_entity_id(EntityHandle this_ent) = 0;
  virtual unsigned int num_entities(int dimension) = 0;

  // ***************************************************************************
  // SECTION IV: Handling DagMC settings
  // ***************************************************************************

  /** retrieve faceting tolerance */
  virtual double faceting_tolerance() = 0;

  /** retrieve overlap thickness */
  double overlap_thickness();

  /** retrieve numerical precision */
  double numerical_precision();

  /** Attempt to set a new overlap thickness tolerance, first checking for sanity */
  void set_overlap_thickness(double new_overlap_thickness);

  /** Attempt to set a new numerical precision , first checking for sanity
   *  Use of this function is discouraged; see top of DagMC.cpp
   */
  void set_numerical_precision(double new_precision);

  // ***************************************************************************
  // SECTION V: Metadata handling
  // ***************************************************************************

  virtual ErrorCode parse_properties(const std::vector<std::string>& keywords,
                                     const std::map<std::string,
                                     std::string>& synonyms = no_synonyms,
                                     const char* delimiters = "_") = 0;
  virtual bool has_prop(EntityHandle eh, const std::string& prop) = 0;
  virtual ErrorCode prop_value(EntityHandle eh, const std::string& prop, std::string& value) = 0;
  virtual ErrorCode detect_available_props(std::vector<std::string>& keywords_out, const char* delimiters = "_") = 0;
  virtual ErrorCode prop_values(EntityHandle eh, const std::string& prop,
                                std::vector<std::string>& value) = 0;

  // ***************************************************************************
  // SECTION VI: Other
  // ***************************************************************************

  /** Write to file */
  ErrorCode write_mesh(const char* ffile);

  virtual ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]) = 0;
  virtual ErrorCode getobb(EntityHandle volume, double center[3], double axis1[3],
                           double axis2[3], double axis3[3]) = 0;
  virtual ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle& root) = 0;

 protected:

  // ***************************************************************************
  // Protected member data
  // ***************************************************************************

  // Generic abstract interface to the mesh
  std::shared_ptr<MeshInterfaceBase> mesh_interface;

  // Generic object to do ray tracing
  std::unique_ptr<RayTracer> ray_tracer;

  // CPointer to class to facilitate error handling
  std::unique_ptr<ErrorHandler> errHandler;

  /** empty synonym map to provide as a default argument to parse_properties() */
  static const std::map<std::string, std::string> no_synonyms;

};

}; // namespace DAGMC

#endif
