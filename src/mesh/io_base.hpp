#ifndef DAG_IOBASE_HPP
#define DAG_IOBASE_HPP

#include "mesh_interface.hpp"

namespace DAGMC {

class IOBase {

 public:

  IOBase() : isOpen(false) {};

  ~IOBase() {
    if (isOpen) {
      close();
    }
  };

  enum Mode { READ, WRITE };

  // Open file for read/write
  virtual bool open(std::string filename, Mode mode) = 0;

  // Read and set mesh attributes. Return success or failure.
  virtual bool readAttributes(MeshAttributes& attributes) = 0;

  // Read and set mesh data.
  virtual bool read() = 0;

  // Write mesh data to file.
  virtual bool write() = 0;

  // Close file
  virtual void close() {};

 protected:

  // Record whether the file is open
  bool isOpen;

};




}




#endif
