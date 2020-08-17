#include <memory>
#include <vector>

namespace DAGMC {

// Abstract base class
class TreeNode : public  std::enable_shared_from_this<TreeNode> {

 public:

  TreeNode(std::shared_ptr<TreeNode> parentIn = nullptr) :
    parent(parentIn) {};
  ~TreeNode() {};

  virtual bool isRoot() {return (parent == nullptr);}
  virtual bool isLeaf() {return (children.empty());}

  void split() {
    // Try to create children from this node
    // Keep going until we reach stop condition
    if (!setChildren() || isLeaf())
      return;
    for (auto& child : children) {
      child->split();
    }
  };

 protected:

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

  // Build the tree
  void build() {
    if (root != nullptr)
      root->split();
  };

 protected:

  // Pointer to the root of the tree
  std::shared_ptr<TreeNode> root;

};


}
