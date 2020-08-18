#ifndef DAG_CONTAINER_HPP
#define DAG_CONTAINER_HPP

#include<set>

// Libmesh headers
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"

// The point of these classes is to abstract away the data types of how we
// actually refer a set of mesh elements, so that they can be iterated
// while being agnostic to  the actual storage mechanism

namespace DAGMC {

typedef libMesh::MeshBase::const_element_iterator const_element_iterator;

// An  abstract external iterator for a Container object.
// Intended for local use only
// In principle have multiple iterators pointing
// at same container
template < class T >
class LocalIterator {

 public:
  LocalIterator() : it_is_void(true) {};
  ~LocalIterator() {};

  // Send internal iterator to start
  virtual void reset() {
    it_is_void = true;
  }

  // Fetch pointers to first and last in list
  virtual T* first() = 0;
  virtual T* last() = 0;

  // Set input pointer to next in list, incrementing internal iterator
  virtual bool getNext(T*& refptr) = 0;

 protected:
  bool it_is_void;
};

// An iterator for libMesh elements
typedef LocalIterator<const libMesh::Elem> ElemIterator;

// Iterator using Libmesh type
class ElemLMIterator : public ElemIterator {

 public:

  ElemLMIterator(const_element_iterator elBegin,
                 const_element_iterator elEnd) :
    elBeginLocal(elBegin),
    elEndLocal(elEnd),
    elIt(elEnd) {};

  ~ElemLMIterator() {};

  const libMesh::Elem* first() override { return &** elBeginLocal; };
  const libMesh::Elem* last() override { return &** elEndLocal; };

  bool getNext(const libMesh::Elem*& next) override {
    if (it_is_void) {
      elIt = elBeginLocal;
      it_is_void = false;
    } else
      elIt++;
    if (elIt != elEndLocal) {
      next = &** elIt;
      return true;
    } else
      return false;
  };

 private:
  // Copies of the beginning and end
  const_element_iterator elBeginLocal;
  const_element_iterator elEndLocal;
  // The actual iterator
  const_element_iterator elIt;

};

//   Iterator over a set of elem pointers
class ElemSetIterator : public ElemIterator {

 public:

  typedef std::set< const libMesh::Elem*>::const_iterator const_set_iterator;

  ElemSetIterator(const_set_iterator elBegin,
                  const_set_iterator elEnd) :
    elBegLocal(elBegin),
    elEndLocal(elEnd) {};

  const libMesh::Elem* first() override { return *elBegLocal; };
  const libMesh::Elem* last() override { return *elEndLocal; };

  bool getNext(const libMesh::Elem*& next) override {
    if (it_is_void) {
      it = elBegLocal;
      it_is_void = false;
    } else
      it++;
    if (it != elEndLocal) {
      next = *it;
      return true;
    } else
      return false;
  }

 private:

  // Copies of the beginning and end
  const_set_iterator elBegLocal;
  const_set_iterator elEndLocal;
  // Iterator
  const_set_iterator it;

};

// An abstract iterable set of type T objects
template < class T >
class Container {

 public:
  Container() {};
  ~Container() {};

  // Create and return an iterator for this container
  // -> clever pointer will deallocate memory once it goes out of scope
  virtual std::shared_ptr<LocalIterator<T> > getIterator() const = 0;

  // // Set input pointer to next in list, incrementing internal iterator
  // virtual bool getNext(T*& refptr) = 0;

};

// A container for libMesh elements
typedef Container<const libMesh::Elem> ElemContainer;

// A container for libMesh elements  constructed from iterators
class ElemConstItContainer : public ElemContainer {

 public:
  ElemConstItContainer(const_element_iterator elemBegin,
                       const_element_iterator elemEnd) :
    elBegin(elemBegin),
    elEnd(elemEnd)
  {};

  ~ElemConstItContainer() {};

  std::shared_ptr<ElemIterator> getIterator() const override {
    std::shared_ptr<ElemIterator> it
      = std::make_shared<ElemLMIterator>(elBegin, elEnd);
    return it;
  };

 private:

  const_element_iterator elBegin;
  const_element_iterator elEnd;

};

// A container for libMesh elements  constructed from a set of pointers
class ElemConstPtrContainer : public ElemContainer {

 public:
  ElemConstPtrContainer(std::set< const libMesh::Elem*> elemsIn) :
    elems(elemsIn) {};

  std::shared_ptr<ElemIterator> getIterator() const override {
    std::shared_ptr<ElemIterator> it
      = std::make_shared<ElemSetIterator>(elems.begin(),
                                          elems.end());
    return it;
  };

 private:

  // Set of constant pointers to contant elements
  const std::set< const libMesh::Elem*> elems;

};

}

#endif
