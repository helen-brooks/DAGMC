//This file was created by H Brooks on 26/06/20
#ifndef DAGMCBASE_HPP

namespace DAGMC {

// This is a base class for the DAGMC object whose derivatives will be MOAB- or
// LibMesh- dependent respectively.
class DagMCBase {
 public:
  // Constructor
  DagMCBase();
  // Destructor
  ~DagMCBase();

  //@TODO Add some pure virtual methods here...
}

} // end DAGMC namespace

#define DAGMCBASE_HPP
#endif