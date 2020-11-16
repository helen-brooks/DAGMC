#ifndef DAG_EXODUSIO_HPP
#define DAG_EXODUSIO_HPP

#ifdef LIBMESH

#include "IO_Base.hpp"
#include "libmesh/exodusII.h"

namespace DAGMC {

class IOBase;

// Currently this is only used for reading attributes.
// But if we later want to in the mesh, could pass a pointer
// to MeshType
class ExodusAttributeReader : public IOBase {

 public:

  ExodusAttributeReader() :
    ex_version(0.),
    io_ws(0),
    comp_ws(8),
    ex_id(0),
    isInit(false) {};
  ~ExodusAttributeReader() {};

  virtual bool open(std::string filename, IOBase::Mode mode);
  virtual bool readAttributes(MeshAttributes& attributes);
  virtual void close();

  // Dummy reader / writer
  virtual bool read() { return false; };
  virtual bool write() { return false; };

 private:

  // Private methods
  bool init();

  // Private data members

  // Word size in bytes
  int io_ws;
  int comp_ws;

  // Version of Exodus: to  be set on reading
  float ex_version;

  // ID code for open file
  int ex_id;

  // Place to store initialisation params
  bool isInit;
  ex_init_params params;


};
}

#endif

#endif
