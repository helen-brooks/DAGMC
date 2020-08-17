#ifndef DAG_CONTAINER_HPP
#define DAG_CONTAINER_HPP

#include<set>

// Libmesh headers
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"

namespace DAGMC {

typedef libMesh::MeshBase::const_element_iterator const_element_iterator;

// An abstract iterable set of type T objects
template < class T >
class Container {

 public:
  Container() {};
  ~Container() {};

  // Send internal iterator to start
  virtual void reset() = 0;

  // Set input pointer to next in list, incrementing internal iterator
  virtual bool getNext(T*& refptr) = 0;

  // Fetch copies of first and last in list
  virtual T* first() = 0;
  virtual T* last() = 0;

};

// A container for libMesh elements
typedef Container<const libMesh::Elem> ElemContainer;

// A container for libMesh elements  constructed from iterators
class ElemConstItContainer : public ElemContainer {

 public:
  ElemConstItContainer(const_element_iterator elemBegin,
                       const_element_iterator elemEnd) :
    elBegin(elemBegin),
    elEnd(elemEnd),
    elIt(elemEnd),
    it_is_void(true) {};

  ~ElemConstItContainer() {};

  void reset() override {
    it_is_void = true;
  }

  const libMesh::Elem* first() override { return &** elBegin; };
  const libMesh::Elem* last() override { return &** elEnd; };

  bool getNext(const libMesh::Elem*& next) override {
    if (it_is_void) {
      elIt = elBegin;
      it_is_void = false;
    } else
      elIt++;
    if (elIt != elEnd) {
      next = &** elIt;
      return true;
    } else
      return false;
  }

 private:

  const_element_iterator elBegin;
  const_element_iterator elEnd;
  const_element_iterator elIt;
  bool it_is_void;

};

// A container for libMesh elements  constructed from a set of pointers
class ElemConstPtrContainer : public ElemContainer {

 public:
  ElemConstPtrContainer(std::set< const libMesh::Elem*> elemsIn) :
    elems(elemsIn) {
    reset();
  };

  void reset() override {
    it_is_void = true;
  }

  const libMesh::Elem* first() override { return *elems.begin(); };
  const libMesh::Elem* last() override { return *elems.end(); };

  bool getNext(const libMesh::Elem*& next) override {
    if (it_is_void) {
      it = elems.begin();
      it_is_void = false;
    } else
      it++;
    if (it != elems.end()) {
      next = *it;
      return true;
    } else
      return false;
  }


 private:

  // Set of constant pointers to contant elements
  const std::set< const libMesh::Elem*> elems;

  // Iterator
  std::set< const libMesh::Elem*>::const_iterator it;
  bool it_is_void;

};

}

#endif
