#ifndef DAG_TREE_HPP
#define DAG_TREE_HPP

#include <memory>
#include <vector>

namespace DAGMC {

// Abstract base class
class TreeNode : public  std::enable_shared_from_this<TreeNode> {

  friend class Tree;

 public:

  TreeNode(std::shared_ptr<TreeNode> parentIn = nullptr) :
    parent(parentIn) {};
  ~TreeNode() {};

  bool isRoot() {return (parent == nullptr);}
  bool isLeaf() {return (children.empty());}

  // Return information on success of construction
  virtual bool isConstructed() const = 0;

  // Return a copy of the children of this node
  std::vector<std::shared_ptr<TreeNode> > getChildren() {
    return children;
  }

 protected:

  // Split is only intended to be called from Tree
  void split() {
    // Try to create children from this node
    // Keep going until we reach stop condition
    if (!setChildren() || isLeaf())
      return;
    for (auto& child : children) {
      child->split();
    }
  };

  // This method defines how the node with be split.
  // Must define a stop condition
  virtual bool setChildren() = 0;

  // Pointers to parent and children nodes. May be null.
  std::shared_ptr<TreeNode> parent;
  std::vector<std::shared_ptr<TreeNode> > children;
};

// Generic Tree that recursively builds itself
class Tree {

 public:

  Tree() {};
  ~Tree() {};

  // Return a copy of the root node
  std::shared_ptr<TreeNode> getRoot() { return root;};

 protected:

  // Recursively build the tree:
  // - intended to be called in derived constructor
  void build() {
    if (root != nullptr)
      root->split();
  };

  // Pointer to the root of the tree
  std::shared_ptr<TreeNode> root;

};

}
#endif
