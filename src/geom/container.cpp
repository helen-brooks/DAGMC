#include "container.hpp"
#include <typeinfo>

namespace DAGMC {


//Helper function
template< class T1, class T2 >
bool compareTypes(T1& t1, T2& t2) {

  if (typeid(t1) == typeid(t2))
    return true;
  else
    return false;

};


bool ElemConstItContainer::checkIfValid() const {

  // Get the underlying data members
  const IterBase* beg_data    = elBegin.data;
  const IterBase* beg_enddata = elBegin.end;
  const IterBase* end_data    = elEnd.data;
  const IterBase* end_enddata = elEnd.end;
  const PredBase* beg_pred    = elBegin.pred;
  const PredBase* end_pred    = elEnd.pred;

  // Check that nothing is null
  if (beg_data == nullptr
      || beg_enddata == nullptr
      || beg_pred == nullptr
      || end_data == nullptr
      || end_enddata == nullptr
      || end_pred == nullptr)
    return false;

  // Check that end points same
  if (!(beg_enddata->equal(end_enddata))) {
    return false;
  }

  // Check unordered
  if (beg_data->equal(end_enddata) &&
      !end_data->equal(end_enddata)) {
    return false;
  }


  // Check unordred for non-trivial endpoints
  if (!beg_data->equal(end_enddata) &&
      !end_data->equal(end_enddata)) {
    libMesh::dof_id_type idBeg = (**beg_data)->id();
    libMesh::dof_id_type idEnd = (**end_data)->id();
    if (idBeg > idEnd) {
      return false;
    }
  }

  // Check that predicates pointed to are the same type
  if (!compareTypes(*beg_pred, *end_pred)) {
    return false;
  }

  // Check that each other's predicates are true
  if (!end_data->equal(end_enddata)) {
    if (!(*beg_pred)(end_data)) {
      return false;
    }
  }
  if (!(*end_pred)(beg_data)) {
    return false;
  }

  // Got to here: passed all checks
  return true;
}


}
