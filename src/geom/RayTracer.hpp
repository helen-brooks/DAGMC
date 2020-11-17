#ifndef DAG_RAY_TRACER
#define DAG_RAY_TRACER

#include "Common.hpp"
#include "Types.hpp"

namespace DAGMC {

class RayTracer {

 public:

  RayTracer() {};
  ~RayTracer() {};

  virtual ErrorCode init()  = 0;

  virtual ErrorCode init_implicit_complement() = 0;

  virtual ErrorCode init_obb()  = 0;

  virtual bool is_implicit_complement(EntityHandle volume)  = 0;

  virtual ErrorCode next_vol(EntityHandle surface, EntityHandle old_volume,
                             EntityHandle& new_volume)  = 0;


  virtual ErrorCode point_in_volume(const EntityHandle volume,
                                    const double xyz[3],
                                    int& result,
                                    const double* uvw,
                                    const moab::GeomQueryTool::RayHistory* history) = 0;

  virtual ErrorCode point_in_volume_slow(EntityHandle volume,
                                         const double xyz[3],
                                         int& result) = 0;

  virtual ErrorCode test_volume_boundary(const EntityHandle volume,
                                         const EntityHandle surface,
                                         const double xyz[3],
                                         const double uvw[3],
                                         int& result,
                                         const moab::GeomQueryTool::RayHistory* history = 0) = 0;

  virtual ErrorCode ray_fire(const EntityHandle volume,
                             const double point[3],
                             const double dir[3],
                             EntityHandle& next_surf,
                             double& next_surf_dist,
                             moab::GeomQueryTool::RayHistory* history = 0,
                             double user_dist_limit = 0,
                             int ray_orientation = 1,
                             moab::OrientedBoxTreeTool::TrvStats* stats = NULL) = 0;

  virtual ErrorCode  closest_to_location(EntityHandle volume,
                                         const double point[3],
                                         double& result,
                                         EntityHandle* closest_surf = 0) = 0;

  virtual ErrorCode get_normal(EntityHandle surf,
                               const double loc[3],
                               double angle[3],
                               const moab::GeomQueryTool::RayHistory* history = 0) = 0;

  virtual ErrorCode  measure_volume(EntityHandle volume,
                                    double& result) = 0;

  virtual ErrorCode  measure_area(EntityHandle surface,
                                  double& result) = 0;

  virtual ErrorCode surface_sense(EntityHandle volume, int num_surfaces,
                                  const EntityHandle* surfaces, int* senses_out) = 0;

  virtual ErrorCode surface_sense(EntityHandle volume, EntityHandle surface,
                                  int& sense_out) = 0;

  virtual double get_numerical_precision() = 0;

  virtual double get_overlap_thickness() = 0;

  virtual void set_numerical_precision(double val) = 0;

  virtual void set_overlap_thickness(double val) = 0;

};
}

#endif
