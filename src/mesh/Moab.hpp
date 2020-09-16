#ifndef DAG_MOAB_HPP
#define DAG_MOAB_HPP

//MOAB headers
#include "MBTagConventions.hpp"
#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/FileOptions.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/GeomQueryTool.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

namespace DAGMC {
// Some typedefs for commonly used moab types
typedef moab::Core Core;
typedef moab::Interface Interface;
typedef moab::TagType TagType;
typedef moab::DataType DataType;
typedef moab::GeomTopoTool GeomTopoTool;
typedef moab::GeomQueryTool GeomQueryTool;
typedef moab::CartVect CartVect;
typedef moab::EntityHandle EntityHandle;
typedef moab::GeomQueryTool::RayHistory RayHistory;
typedef moab::OrientedBoxTreeTool OrientedBoxTreeTool;
typedef moab::Range Range;
typedef moab::Tag Tag;
}

#endif





