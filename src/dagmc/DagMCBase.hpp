//This file was created by H Brooks on 26/06/20
#ifndef DAGMCBASE_HPP
#define DAGMCBASE_HPP

// DAGMC headers
#include "DagMCVersion.hpp"
#include "common.hpp"
#include "Moab.hpp"
#include "Error.hpp"

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

  virtual ErrorCode load_file(const char* cfile) = 0;
  virtual ErrorCode load_existing_contents() = 0;
  virtual ErrorCode init_OBBTree() = 0;
  virtual ErrorCode setup_impl_compl() = 0;
  virtual ErrorCode setup_obbs() = 0;
  virtual ErrorCode setup_indices() = 0;

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
  virtual ErrorCode detect_available_props(std::vector<std::string>& keywords_out, const char* delimiters = "_") = 0;
  virtual ErrorCode prop_values(EntityHandle eh, const std::string& prop,
                                std::vector<std::string>& value) = 0;

  // ***************************************************************************
  // SECTION VI: Other
  // ***************************************************************************

  virtual ErrorCode write_mesh(const char* ffile) = 0;
  virtual ErrorCode getobb(EntityHandle volume, double minPt[3], double maxPt[3]) = 0;
  virtual ErrorCode getobb(EntityHandle volume, double center[3], double axis1[3],
                           double axis2[3], double axis3[3]) = 0;
  virtual ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle& root) = 0;

 protected:

  // ***************************************************************************
  // Protected member data
  // ***************************************************************************

  std::unique_ptr<ErrorHandler> errHandler;
  double facetingTolerance;
  /** empty synonym map to provide as a default argument to parse_properties() */
  static const std::map<std::string, std::string> no_synonyms;

};

}; // namespace DAGMC

#endif
