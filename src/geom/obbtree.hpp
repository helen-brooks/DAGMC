#ifndef DAG_OBBTREE_HPP
#define DAG_OBBTREE_HPP

#include "tree.hpp"
#include "obb.hpp"

namespace DAGMC {

class OBBTree : public Tree {

 public:

  OBBTree(const_element_iterator elemBegin,
          const_element_iterator elemEnd) {

    // Create a Bounding Box as the root node.
    root = std::make_shared<OrientedBoundingBox>(elemBegin, elemEnd);

    // Recursively build the tree.
    build();

  };

  ~OBBTree() {}
};
}

#endif
