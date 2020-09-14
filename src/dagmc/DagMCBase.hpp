//This file was created by H Brooks on 26/06/20
#ifndef DAGMCBASE_HPP
#define DAGMCBASE_HPP

// Header files

// DAGMC headers
#include "DagMCVersion.hpp"

//MOAB headers
// TODO remove moab headers
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/GeomQueryTool.hpp"

// Standard library headers
#include <assert.h>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
//#include <algorithm>
#include <set>
#include <climits>

//#include <ctype.h>
//#include <string.h>
//#include <stdlib.h>
//#include <stdio.h>

//#include <math.h>
//#ifndef M_PI  /* windows */
//# define M_PI 3.14159265358979323846
//#endif

namespace DAGMC {

// Backwards compatible error code handling
#ifdef MOAB_TYPES_HPP
enum ErrorCode {
  DAG_SUCCESS = moab::MB_SUCCESS,
  DAG_INDEX_OUT_OF_RANGE = moab::MB_INDEX_OUT_OF_RANGE,
  DAG_TYPE_OUT_OF_RANGE = moab::MB_TYPE_OUT_OF_RANGE,
  DAG_MEMORY_ALLOCATION_FAILED = moab::MB_MEMORY_ALLOCATION_FAILED,
  DAG_ENTITY_NOT_FOUND = moab::MB_ENTITY_NOT_FOUND,
  DAG_MULTIPLE_ENTITIES_FOUND = moab::MB_MULTIPLE_ENTITIES_FOUND,
  DAG_TAG_NOT_FOUND = moab::MB_TAG_NOT_FOUND,
  DAG_FILE_DOES_NOT_EXIST = moab::MB_FILE_DOES_NOT_EXIST,
  DAG_FILE_WRITE_ERROR = moab::MB_FILE_WRITE_ERROR,
  DAG_NOT_IMPLEMENTED = moab::MB_NOT_IMPLEMENTED,
  DAG_ALREADY_ALLOCATED = moab::MB_ALREADY_ALLOCATED,
  DAG_VARIABLE_DATA_LENGTH = moab::MB_VARIABLE_DATA_LENGTH,
  DAG_INVALID_SIZE = moab::MB_INVALID_SIZE,
  DAG_UNSUPPORTED_OPERATION = moab::MB_UNSUPPORTED_OPERATION,
  DAG_UNHANDLED_OPTION = moab::MB_UNHANDLED_OPTION,
  DAG_STRUCTURED_MESH = moab::MB_STRUCTURED_MESH,
  DAG_FAILURE = moab::MB_FAILURE
};
#else
//Hard code to take same implicit values as above
enum ErrorCode {
  DAG_SUCCESS = 0,
  DAG_INDEX_OUT_OF_RANGE,
  DAG_TYPE_OUT_OF_RANGE,
  DAG_MEMORY_ALLOCATION_FAILED,
  DAG_ENTITY_NOT_FOUND,
  DAG_MULTIPLE_ENTITIES_FOUND,
  DAG_TAG_NOT_FOUND,
  DAG_FILE_DOES_NOT_EXIST,
  DAG_FILE_WRITE_ERROR,
  DAG_NOT_IMPLEMENTED,
  DAG_ALREADY_ALLOCATED,
  DAG_VARIABLE_DATA_LENGTH,
  DAG_INVALID_SIZE,
  DAG_UNSUPPORTED_OPERATION,
  DAG_UNHANDLED_OPTION,
  DAG_STRUCTURED_MESH,
  DAG_FAILURE
};
#endif

// Some typedefs for commonly used moab types
// @TODO Generalise and reduce MOAB dependence
typedef moab::EntityHandle EntityHandle;
typedef moab::GeomQueryTool::RayHistory RayHistory;
typedef moab::OrientedBoxTreeTool OrientedBoxTreeTool;
typedef moab::Range Range;
typedef moab::Tag Tag;

static const int vertex_handle_idx = 0;
static const int curve_handle_idx = 1;
static const int surfs_handle_idx = 2;
static const int vols_handle_idx = 3;
static const int groups_handle_idx = 4;

/** Helper class to handle error codes */
class ErrorHandler {
 public:
  ErrorHandler() : _code(DAG_SUCCESS) {};
  ~ErrorHandler() {};
  /** Check the error code for a set error.
   *  Print message if error is detected. */
  virtual void checkSetErr(ErrorCode rval, std::string msg) = 0;

#ifdef MOAB_TYPES_HPP
  /** Overloaded version to cast moab ErrorCode type as DAGMC native ErrorCode */
  void checkSetErr(moab::ErrorCode mbcode, std::string msg) {
    checkSetErr(ErrorCode(mbcode), msg);
  };
#endif
  /** Retrieve the latest error code */
  ErrorCode code() { return _code; };

 protected:
  ErrorCode _code;
};

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
 *   2) DAG->setup_obb();
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

  /** Get subversion revision of this file (DagMC.hpp) */
  static unsigned int interface_revision();

  // ***************************************************************************
  // SECTION I: Geometry Initialization
  // ***************************************************************************

  virtual ErrorCode load_file(const char* cfile) = 0;
  virtual ErrorCode setup_impl_compl() = 0;
  /** \brief initializes the geometry and OBB tree structure for ray firing acceleration
  *
  * This method can be called after load_file to fully initialize DAGMC.  It calls
  * methods to set up the geometry, create the implicit complement, generate an
  * OBB tree from the faceted representation of the geometry, and build the
  * cross-referencing indices.
  */
  ErrorCode init_OBBTree();
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
  /** Thin wrapper around build_indices().*/
  ErrorCode setup_indices();

  // ***************************************************************************
  // SECTION II: Fundamental Geometry Operations/Queries
  // ***************************************************************************

  virtual ErrorCode ray_fire(const EntityHandle volume, const double ray_start[3],
                             const double ray_dir[3], EntityHandle& next_surf,
                             double& next_surf_dist,
                             RayHistory* history = NULL,
                             double dist_limit = 0, int ray_orientation = 1,
                             OrientedBoxTreeTool::TrvStats* stats = NULL) = 0;
  virtual ErrorCode point_in_volume(const EntityHandle volume, const double xyz[3],
                                    int& result, const double* uvw = NULL,
                                    const RayHistory* history = NULL) = 0;
  virtual ErrorCode point_in_volume_slow(const EntityHandle volume, const double xyz[3],
                                         int& result) = 0;
  virtual ErrorCode test_volume_boundary(const EntityHandle volume,
                                         const EntityHandle surface,
                                         const double xyz[3], const double uvw[3],
                                         int& result,
                                         const RayHistory* history = NULL) = 0;
  virtual ErrorCode closest_to_location(EntityHandle volume, const double point[3],
                                        double& result, EntityHandle* surface = 0) = 0;
  virtual ErrorCode measure_volume(EntityHandle volume, double& result) = 0;
  virtual ErrorCode measure_area(EntityHandle surface, double& result) = 0;
  virtual ErrorCode surface_sense(EntityHandle volume, int num_surfaces,
                                  const EntityHandle* surfaces, int* senses_out) = 0;
  virtual ErrorCode surface_sense(EntityHandle volume, EntityHandle surface,
                                  int& sense_out) = 0;
  virtual ErrorCode get_angle(EntityHandle surf, const double xyz[3], double angle[3],
                              const RayHistory* history = NULL) = 0;
  virtual ErrorCode next_vol(EntityHandle surface, EntityHandle old_volume,
                             EntityHandle& new_volume) = 0;

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

  virtual double overlap_thickness() = 0;
  virtual void set_overlap_thickness(double new_overlap_thickness) = 0;
  /** retrieve faceting tolerance */
  double faceting_tolerance() { return facetingTolerance; }

  // ***************************************************************************
  // SECTION V: Metadata handling
  // ***************************************************************************

  virtual ErrorCode parse_properties(const std::vector<std::string>& keywords,
                                     const std::map<std::string,
                                     std::string>& synonyms = no_synonyms,
                                     const char* delimiters = "_") = 0;
  virtual bool has_prop(EntityHandle eh, const std::string& prop) = 0;
  virtual bool is_implicit_complement(EntityHandle volume) = 0;
  virtual ErrorCode prop_value(EntityHandle eh, const std::string& prop, std::string& value) = 0;

  /** Detect all the property keywords that appear in the loaded geometry
   *
   *  @param keywords_out The result list of keywords.  This list could be
   *        validly passed to parse_properties().
   */
  ErrorCode detect_available_props(std::vector<std::string>& keywords_out, const char* delimiters = "_");

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
                        std::vector<std::string>& value);

  // SECTION VI: Other
  virtual ErrorCode write_mesh(const char* ffile, const int flen) = 0;
  virtual ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]) = 0;
  virtual ErrorCode getobb(EntityHandle volume, double center[3], double axis1[3],
                           double axis2[3], double axis3[3]) = 0;
  virtual ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle& root) = 0;

 protected:
  // Typedefs
  /** a common type within the property and group name functions */
  typedef std::map<std::string, std::string> prop_map;

  // ***************************************************************************
  // Protected Methods
  // ***************************************************************************

  /** \brief Loading code shared by load_file and load_existing_contents
   *  Called by load_file and load_existing_contents
  */
  virtual ErrorCode finish_loading() = 0;

  virtual ErrorCode setup_obbs() = 0;

  virtual ErrorCode find_geomsets() = 0;

  /** \brief sets up ranges of the volume and surface entity sets
  *
  * Helper function for setup_indices. Sets ranges containing
  * all volumes and surfaces.Called by setup_indices
  */
  virtual ErrorCode setup_geometry(Range& surfs, Range& vols) = 0;
  /** \brief build internal index vectors that speed up handle-by-id, etc.
  * Called by setup_indices
  */
  virtual ErrorCode build_indices(Range& surfs, Range& vols) = 0;
  /** \brief Convert a property tag's value on a handle to a list of strings
  * Called by prop_values
  */
  virtual ErrorCode unpack_packed_string(Tag tag, EntityHandle eh,
                                         std::vector<std::string>& values) = 0;
  /** \brief Store the name of a group in a string
  * Called by parse_group_name
  */
  virtual ErrorCode get_group_name(EntityHandle group_set, std::string& name) = 0;

  /** \brief tokenize the metadata stored in group names
  * - basically borrowed from ReadCGM.cpp.
  * Called by parse_group_name
  * */
  void tokenize(const std::string& str, std::vector<std::string>& tokens,
                const char* delimiters = "_") const;

  /** \brief Parse a group name into a set of key:value pairs
  * Called by detect_available_props, parse_properties
  */
  ErrorCode parse_group_name(EntityHandle group_set, prop_map& result, const char* delimiters = "_");

  std::vector<EntityHandle>& surf_handles() {
    return entHandles[surfs_handle_idx];
  };
  std::vector<EntityHandle>& vol_handles() {
    return entHandles[vols_handle_idx];
  };
  std::vector<EntityHandle>& group_handles() {
    return entHandles[groups_handle_idx];
  };

  // ***************************************************************************
  // Protected member data
  // ***************************************************************************
  /** Pointer to ErrorHandler object */
  std::unique_ptr<ErrorHandler> errHandler;
  /** empty synonym map to provide as a default argument to parse_properties() */
  static const std::map<std::string, std::string> no_synonyms;
  /** map from the canonical property names to the tags representing them */
  std::map<std::string, Tag> property_tagmap;
  /** store some lists indexed by handle */
  std::vector<EntityHandle> entHandles[5];
  double facetingTolerance;
};

}; // namespace DAGMC

#endif
