#ifndef DAG_IO_HPP
#define DAG_IO_HPP

#ifdef LIBMESH

#include "ExodusII_IO.hpp"

namespace DAGMC {

//  IO factory (in case we add more file types later)
std::shared_ptr<IOBase> getIOPtr(std::string filename) {

  std::shared_ptr<IOBase> ioptr(nullptr);
  if (filename.rfind(".e") ==  std::string::npos) {
    std::cerr << "Unrecognised mesh file type with name "
              << filename << std::endl;
  } else {
    ioptr  = std::make_shared<ExodusAttributeReader>();
  }
  return ioptr;

}

}

#endif

#endif
