#include "exodusII_io.hpp"

#ifdef LIBMESH
namespace DAGMC {

bool ExodusAttributeReader::open(std::string filename, IOBase::Mode mode) {

  // If we already have a file open don't try to open another
  // User should close first
  if (isOpen)
    return false;

  int EX_MODE = (mode == IOBase::READ) ? EX_READ : EX_WRITE;
  ex_id = ex_open(filename.c_str(), EX_MODE, &comp_ws, &io_ws, &ex_version);
  if (ex_id < 0)
    return false;
  else {
    isOpen = true;
    return true;
  }
}

bool ExodusAttributeReader::readAttributes(MeshAttributes& attributes) {

  // If nothing is open cannot read
  if (!isOpen)
    return false;

  if (!init())
    return false;

  std::vector<ex_entity_id> block_ids;
  block_ids.resize(params.num_elem_blk);

  int ex_ret = ex_get_ids(ex_id, EX_ELEM_BLOCK, block_ids.data());
  if (ex_ret < 0)
    return false;

  for (auto& block : block_ids) {

    int num_elems;
    int num_nodes;
    int num_edges;
    int num_faces;
    int num_attrs;
    ex_ret = ex_get_block(ex_id, EX_ELEM_BLOCK, block, 0, &num_elems, &num_nodes,
                          &num_edges, &num_faces, &num_attrs);

    // Expect at least two attributes
    if (ex_ret < 0 || num_attrs < 2)
      return false;

    std::vector<double> block_attr;
    block_attr.resize(num_attrs * num_elems);
    ex_ret =  ex_get_attr(ex_id, EX_ELEM_BLOCK, block,
                          block_attr.data());

    if (ex_ret < 0)
      return false;

    // Get first two attributes: expect to be surface senses
    SurfaceSenses senses;

    // Downcast attributes as ints
    senses.forwards = intmax_t(block_attr.at(0));
    senses.backwards = intmax_t(block_attr.at(1));

    // Block id  = surf id
    intmax_t surf(block);

    // Attempt to save, but don't overwrite existing entries
    if (attributes.senseData.find(surf) != attributes.senseData.end()) {
      return false;
    } else {
      attributes.senseData[surf] = senses;
    }

  }

  // Sanity check
  if (attributes.senseData.size() != params.num_elem_blk) {
    return false;
  }

  // Got to here: success
  return true;

}

void ExodusAttributeReader::close() {

  // If nothing is open nothing to do.
  if (!isOpen)
    return;
  else {
    ex_close(ex_id);
  }
  return;
}

bool ExodusAttributeReader::init() {

  // Don't attempt to initialise if no open file
  if (!isOpen)
    return false;

  // Don't re-initialise
  if (isInit)
    return true;

  // Fetch initialisation data
  int ex_ret = ex_get_init_ext(ex_id, &params);

  // Record success / failure
  if (ex_ret >= 0) {
    isInit = true;
  }

  return isInit;
}

}
#endif
