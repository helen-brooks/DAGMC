#include "GeomTopoTool.hpp"
#include <libmesh/boundary_info.h>

namespace DAGMC {


bool GeomTopoToolLM::setupGeometry() {

  // TODO: do some checks, e.g. dimensionality of the mesh.

  // Get surfaces = element sets = subdomains
  std::set< libMesh::subdomain_id_type > surf_ids;
  meshPtr->subdomain_ids(surf_ids);

  if (surf_ids.size() == 0) {
    std::cerr << "No surfaces were found in this file." << std::endl;
    return false;
  }

  // Fetch a reference to the boundary info
  const libMesh::BoundaryInfo bInfo = meshPtr->get_boundary_info();

  // Get volumes = nodesets (convention)
  const std::set< libMesh::boundary_id_type > vols
    = bInfo.get_node_boundary_ids();

  // Extract mapping of surfaces to volumes￼
  // 1) Initialise
  for (auto vol : vols) {
    vol_to_surfs[int(vol)] = std::vector<int>();
  }
  vol_to_surfs[IMPLICIT_COMPLEMENT] = std::vector<int>();

  // 2) Loop over surfaces
  for (auto isurf : surf_ids) {

    // Get first element in surface
    const libMesh::Elem& elem =
        **(meshPtr->active_local_subdomain_elements_begin(isurf));

    std::set<libMesh::boundary_id_type> shared_vols;

    // Loop over the nodes in element
    unsigned int nNodes = elem.n_nodes();
    for (unsigned int iNode = 0; iNode < nNodes; iNode++) {
      const libMesh::Node* nodePtr = elem.node_ptr(iNode);

      auto shared_vols_so_far = shared_vols;
      shared_vols.clear();

      // Get all boundary IDs for this node
      std::vector< libMesh::boundary_id_type > bIDs;
      bInfo.boundary_ids(nodePtr, bIDs);

      // Save shared volume IDs
      unsigned int nVolsNode = 0;
      for (auto ID : bIDs) {
        // Check this is volume ID
        if (vols.find(ID) != vols.end()) {
          nVolsNode++;
          //If first node, just save vol
          if (iNode == 0) {
            shared_vols.insert(ID);
          } else {
            // Only save if this vol is shared with all others so far.
            if (shared_vols_so_far.find(ID) != shared_vols_so_far.end()) {
              shared_vols.insert(ID);
            }
          }
        }
      }

      // Check each node had at least one vol
      if (nVolsNode < 1) {
        //Error, exit
        return false;
      }
      // Check that # shared nodes is at least one at all times
      else if (shared_vols.size() < 1) {
        //Error, exit
        return false;
      }
    } // End loop over nodes

    // We now have the volumes for this surface.
    // Should be identically two, unless bounded by IC
    if (shared_vols.size() != 2) {
      if (shared_vols.size() == 1) {
        shared_vols.insert(IMPLICIT_COMPLEMENT);
      } else {
        //Error, exit
        return false;
      }
    }

    // Cast surface ID as int
    int surf = int(isurf);

    // Get pair of vols and cast as int
    auto it = shared_vols.begin();
    int vol1 = int(*it);
    it++;
    int vol2 = int(*it);

    if (surf == VOID_INDEX || vol1 == VOID_INDEX || vol2 == VOID_INDEX) {
      std::cerr << "The index " << VOID_INDEX << " is reserved.";
      std::cerr << " Please don't use this as an identifier.";
      std::cerr << std::endl;
      return false;
    }

    //Save
    surf_to_vols[surf] = std::make_pair(vol1, vol2);
    vol_to_surfs[vol1].push_back(surf);
    vol_to_surfs[vol2].push_back(surf);

  }// end loop over surfaces

  // Perhaps optionally check that implicit complement forms a closed convex surface.
  //if(checkIC()) return false;

  // Got to here - success!
  return true;
}; // end setupGeometry

void GeomTopoToolLM::print() {

  std::cout << "Found " << surf_to_vols.size() << " surfaces and "
            << vol_to_surfs.size() << " volumes" << std::endl;
  for (auto& surf2vol : surf_to_vols) {
    std::cout << "Surface : " << surf2vol.first
              << " --> vols (" << surf2vol.second.first
              << " , " << surf2vol.second.second
              << ")" << std::endl;
  }
  for (auto& surfset : vol_to_surfs) {
    std::cout << "Vol : " << surfset.first
              << " --> surfs (";
    for (auto& surf : surfset.second) {
      std::cout << surf << " ";
    }
    std::cout << ")" << std::endl;
  }
}

int GeomTopoToolLM::getVolID(unsigned int index) {

  if (index >= nVols())
    return VOID_INDEX;
  auto volIt = vol_to_surfs.begin();
  std::advance(volIt, index);
  return volIt->first;

}

int GeomTopoToolLM::getSurfID(unsigned int index) {

  if (index >= nSurfs())
    return VOID_INDEX;
  auto surfIt = surf_to_vols.begin();
  std::advance(surfIt, index);
  return surfIt->first;

}

//Get surfaces belonging to volume
const std::vector<int>& GeomTopoToolLM::getSurfs(unsigned int index) {

  int id = getVolID(index);
  // Will throw exception for index out of range
  return vol_to_surfs.at(id);

}

//Get volumes belonging to a surface
const std::pair<int, int>& GeomTopoToolLM::getVolPair(unsigned int index) {

  int id = getSurfID(index);
  // Will throw exception for index out of range
  return surf_to_vols.at(id);

}


}
