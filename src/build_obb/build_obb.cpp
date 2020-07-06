#include <DagMCmoab.hpp>
#include "moab/ProgOptions.hpp"

int main(int argc, char* argv[]) {

  std::string dag_file;
  std::string out_file;
  bool verbose = false;

  ProgOptions po("build_obb: A tool to prebuild your DAGMC OBB Tree");

  po.addOpt<void>("verbose,v", "Verbose output", &verbose);
  po.addRequiredArg<std::string>("dag_file", "Path to DAGMC file to proccess", &dag_file);
  po.addOpt<std::string>("output,o", "Specify the output filename (default "")", &out_file);

  po.addOptionHelpHeading("Options for loading files");

  po.parseCommandLine(argc, argv);

  // make new DagMC
  DAGMC::DagMCBase* DAG = new DAGMC::DagMCmoab();

  DAGMC::ErrorCode rval;

  // sets the output filename if none specified
  if (out_file == "") {
    int pos = dag_file.find(".h5m");
    if (pos != std::string::npos) {
      out_file = dag_file.substr(0, pos);
    } else {
      out_file = dag_file;
    }
    out_file = out_file + "_obb.h5m";
    std::cout << "Setting default outfile to be " << out_file << std::endl;
  }

  // read geometry
  rval = DAG->load_file(dag_file.c_str());
  if (DAGMC::DAG_SUCCESS != rval) {
    std::cerr << "DAGMC failed to read input file: " << dag_file << std::endl;
    exit(EXIT_FAILURE);
  }

  // initialize geometry
  rval = DAG->init_OBBTree();
  if (DAGMC::DAG_SUCCESS != rval) {
    std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  // write the new file
  rval = DAG->write_mesh(out_file.c_str(), out_file.length());
  if (DAGMC::DAG_SUCCESS != rval) {
    std::cerr << "DAGMC failed write file with OBB" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  return 0;
}
