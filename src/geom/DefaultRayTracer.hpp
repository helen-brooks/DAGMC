#ifndef DAG_DEFAULT_RAY_TRACER
#define DAG_DEFAULT_RAY_TRACER

#include "RayTracer.hpp"
#include "MoabInterface.hpp"

#ifdef DOUBLE_DOWN
#include "RTI.hpp"
#include "MOABRay.h"
#endif

namespace DAGMC {
#ifndef DOUBLE_DOWN
using RayTracingInterface = GeomQueryTool;
#endif

class DefaultRayTracer : public RayTracer {

 public:

  DefaultRayTracer(std::shared_ptr<MoabInterface> mesh_interface, double overlap_tolerance, double numerical_precision);
  ~DefaultRayTracer() {};

  ErrorCode init() override;

  ErrorCode init_implicit_complement() override;

  ErrorCode init_obb() override;

  bool is_implicit_complement(EntityHandle volume) override;

  ErrorCode next_vol(EntityHandle surface, EntityHandle old_volume,
                     EntityHandle& new_volume) override;


  ErrorCode point_in_volume(const moab::EntityHandle volume,
                            const double xyz[3],
                            int& result,
                            const double* uvw,
                            const moab::GeomQueryTool::RayHistory* history) override;

  ErrorCode point_in_volume_slow(moab::EntityHandle volume,
                                 const double xyz[3],
                                 int& result) override;

  ErrorCode test_volume_boundary(const moab::EntityHandle volume,
                                 const moab::EntityHandle surface,
                                 const double xyz[3],
                                 const double uvw[3],
                                 int& result,
                                 const moab::GeomQueryTool::RayHistory* history = 0) override;

  ErrorCode ray_fire(const moab::EntityHandle volume,
                     const double point[3],
                     const double dir[3],
                     moab::EntityHandle& next_surf,
                     double& next_surf_dist,
                     moab::GeomQueryTool::RayHistory* history = 0,
                     double user_dist_limit = 0,
                     int ray_orientation = 1,
                     OrientedBoxTreeTool::TrvStats* stats = NULL) override;

  ErrorCode  closest_to_location(moab::EntityHandle volume,
                                 const double point[3],
                                 double& result,
                                 moab::EntityHandle* closest_surf = 0) override;

  ErrorCode get_normal(moab::EntityHandle surf,
                       const double loc[3],
                       double angle[3],
                       const moab::GeomQueryTool::RayHistory* history = 0) override;

  ErrorCode  measure_volume(moab::EntityHandle volume,
                            double& result) override;

  ErrorCode  measure_area(moab::EntityHandle surface,
                          double& result) override;

  ErrorCode surface_sense(EntityHandle volume, int num_surfaces,
                          const EntityHandle* surfaces, int* senses_out) override;

  ErrorCode surface_sense(EntityHandle volume, EntityHandle surface,
                          int& sense_out) override;

  double get_numerical_precision() override;

  double get_overlap_thickness() override;

  void set_numerical_precision(double val) override;

  void set_overlap_thickness(double val) override;

 private:

  std::unique_ptr<ErrorHandler> errHandler;
  std::unique_ptr<RayTracingInterface> RTI;
  std::shared_ptr<moab::GeomTopoTool> GTT;

};
}

#endif
