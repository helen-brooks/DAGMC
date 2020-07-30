// C++ include files that we need
#include <iostream>
// Functions to initialize the library.
#include "libmesh/libmesh.h"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "GeomTopoTool.hpp"

int main (int argc, char ** argv)
{

  // If exodus support is not enabled exit here
#ifndef LIBMESH_HAVE_EXODUS_API
  std::cerr<<"Please re-build libMesh with ExodusII support."<<std::endl;
  return EXIT_FAILURE;
#endif

  if(argc!=2){
    std::cerr<<"Usage: "<< argv[0] <<" [exodus mesh file]"<<std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  if(filename.rfind(".e") >= filename.size()){
    std::cerr<< "Please provide a *.e[xd] (Exodus II) mesh file as input."
	     <<std::endl;      
    std::cerr<<"Usage: "<< argv[0] <<" [exodus mesh file]"<<std::endl;
    return EXIT_FAILURE;
  }
    
  // Initialize the library.
  libMesh::LibMeshInit init (argc, argv);

  std::cout<<"Creating mesh "<<std::endl;

  //Create the mesh 
  std::shared_ptr<libMesh::Mesh> meshPtr
    = std::make_shared<libMesh::Mesh>(init.comm());

  std::cout<<"Loading mesh from "<< filename << std::endl;
  
  // Read the input mesh.
  meshPtr->read(filename);
  
  DAGMC::GeomTopoToolLM gTT(meshPtr);
  if(!gTT.setupGeometry()){
    std::cerr<<"Failed to setup geometry for DAGMC"<<std::endl;
  }
  else{
    std::cout<<"Succeeding in setting up geometry for DAGMC"<<std::endl;
    gTT.print();
  }
  
  return EXIT_SUCCESS;
}

