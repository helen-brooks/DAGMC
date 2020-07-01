//This file was created by H Brooks on 26/06/20
#ifndef DAGMCBASE_HPP
#define DAGMCBASE_HPP

// Header files
// TODO remove moab headers
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/GeomQueryTool.hpp"
#include "DagMCVersion.hpp"

#include <assert.h>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace DAGMC {

// Backwards compatible error code handling
#ifdef MOAB_TYPES_HPP
enum ErrorCode {
  DAG_SUCCESS = moab::MB_SUCCESS,
  DAG_TAG_NOT_FOUND = moab::MB_TAG_NOT_FOUND,
  DAG_FAILURE = moab::MB_FAILURE,
};
#else
//Hard code to take same implicit values as above
enum ErrorCode {
  DAG_SUCCESS = 0,
  DAG_TAG_NOT_FOUND = 6,
  DAG_FAILURE = 16,
};
#endif

// Some typedefs for commonly used moab types
// @TODO Generalise and reduce MOAB dependence
typedef moab::EntityHandle EntityHandle;
typedef moab::GeomQueryTool::RayHistory RayHistory;
typedef moab::Range Range;
typedef moab::OrientedBoxTreeTool::TrvStats TrvStats;
typedef moab::Tag Tag;

static const int vertex_handle_idx = 0;
static const int curve_handle_idx = 1;
static const int surfs_handle_idx = 2;
static const int vols_handle_idx = 3;
static const int groups_handle_idx = 4;

class ErrorHandler {
 public:
  ErrorHandler() : _code(DAG_SUCCESS) {};
  ~ErrorHandler() {};

  virtual void checkSetErr(ErrorCode rval, std::string msg) = 0;

#ifdef MOAB_TYPES_HPP
  //overload functionality if moab exists
  void checkSetErr(moab::ErrorCode mbcode, std::string msg) {
    checkSetErr(ErrorCode(mbcode), msg);
  };
#endif

  ErrorCode code() { return _code; };

 protected:
  ErrorCode _code;
};

// This is a base class for the DAGMC object whose derivatives will be MOAB- or
// LibMesh- dependent respectively.
class DagMCBase {
 public:
  // Constructor
  DagMCBase();
  // Destructor
  ~DagMCBase();

  // Public methods

  // Return the version of this library
  static float version(std::string* version_string = NULL);

  // SECTION I: Geometry Initialization
  virtual ErrorCode load_file(const char* cfile) = 0;
  virtual ErrorCode init_OBBTree() = 0;
  virtual ErrorCode setup_impl_compl() = 0;
  ErrorCode load_existing_contents();
  ErrorCode setup_indices();

  // SECTION II: Fundamental Geometry Operations/Queries
  virtual ErrorCode ray_fire(const EntityHandle volume, const double ray_start[3],
                             const double ray_dir[3], EntityHandle& next_surf,
                             double& next_surf_dist,
                             RayHistory* history = NULL,
                             double dist_limit = 0, int ray_orientation = 1,
                             TrvStats* stats = NULL) = 0;
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

  // SECTION III: Indexing & Cross-referencing
  virtual EntityHandle entity_by_id(int dimension, int id) = 0;
  virtual EntityHandle entity_by_index(int dimension, int index) = 0;
  virtual int id_by_index(int dimension, int index) = 0;
  virtual int index_by_handle(EntityHandle handle) = 0;
  virtual int get_entity_id(EntityHandle this_ent) = 0;
  virtual unsigned int num_entities(int dimension) = 0;

  // SECTION IV: Handling DagMC settings
  virtual double overlap_thickness() = 0;
  virtual void set_overlap_thickness(double new_overlap_thickness) = 0;
  double faceting_tolerance() { return facetingTolerance; }

  // SECTION V: Metadata handling
  virtual ErrorCode parse_properties(const std::vector<std::string>& keywords,
                                     const std::map<std::string,
                                     std::string>& synonyms = no_synonyms,
                                     const char* delimiters = "_") = 0;
  virtual bool has_prop(EntityHandle eh, const std::string& prop) = 0;
  virtual bool is_implicit_complement(EntityHandle volume) = 0;

  ErrorCode detect_available_props(std::vector<std::string>& keywords_out, const char* delimiters = "_");
  ErrorCode prop_values(EntityHandle eh, const std::string& prop,
                        std::vector<std::string>& value);

  // SECTION VI: Other
  virtual ErrorCode write_mesh(const char* ffile, const int flen) = 0;
  virtual ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]) = 0;
  virtual ErrorCode getobb(EntityHandle volume, double center[3], double axis1[3],
                           double axis2[3], double axis3[3]) = 0;

 protected:
  // Typedef
  typedef std::map<std::string, std::string> prop_map;

  // Methods

  // Called by load_file and load_existing_contents
  virtual ErrorCode finish_loading() = 0;

  // Called by setup_indices
  virtual ErrorCode setup_geometry(Range& surfs, Range& vols) = 0;
  virtual ErrorCode build_indices(Range& surfs, Range& vols) = 0;
  // Called by prop_values
  virtual ErrorCode unpack_packed_string(Tag tag, EntityHandle eh,
                                         std::vector<std::string>& values) = 0;
  // Called by parse_group_name
  virtual ErrorCode get_group_name(EntityHandle group_set, std::string& name) = 0;
  void tokenize(const std::string& str, std::vector<std::string>& tokens,
                const char* delimiters = "_") const;

  // Called by detect_available_props, parse_properties
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

  // Members
  std::unique_ptr<ErrorHandler> errHandler;
  double facetingTolerance;
  static const std::map<std::string, std::string> no_synonyms;
  std::vector<EntityHandle> entHandles[5];
  std::map<std::string, Tag> property_tagmap;
};

}; // namespace DAGMC

#endif
