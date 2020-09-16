#ifndef DAG_ERROR_HPP
#define DAG_ERROR_HPP

#include "Moab.hpp"

namespace DAGMC {

// Backwards compatible error code handling
#ifdef MOAB_TYPES_HPP
enum ErrorCode {
  DAG_SUCCESS = moab::MB_SUCCESS,
  DAG_INDEX_OUT_OF_RANGE = moab::MB_INDEX_OUT_OF_RANGE,
  DAG_TYPE_OUT_OF_RANGE = moab::MB_TYPE_OUT_OF_RANGE,
  DAG_MEMORY_ALLOCATION_FAILED = moab::MB_MEMORY_ALLOCATION_FAILED,
  DAG_ENTITY_NOT_FOUND = moab::MB_ENTITY_NOT_FOUND,
  DAG_MULTIPLE_ENTITIES_FOUND = moab::MB_MULTIPLE_ENTITIES_FOUND,
  DAG_TAG_NOT_FOUND = moab::MB_TAG_NOT_FOUND,
  DAG_FILE_DOES_NOT_EXIST = moab::MB_FILE_DOES_NOT_EXIST,
  DAG_FILE_WRITE_ERROR = moab::MB_FILE_WRITE_ERROR,
  DAG_NOT_IMPLEMENTED = moab::MB_NOT_IMPLEMENTED,
  DAG_ALREADY_ALLOCATED = moab::MB_ALREADY_ALLOCATED,
  DAG_VARIABLE_DATA_LENGTH = moab::MB_VARIABLE_DATA_LENGTH,
  DAG_INVALID_SIZE = moab::MB_INVALID_SIZE,
  DAG_UNSUPPORTED_OPERATION = moab::MB_UNSUPPORTED_OPERATION,
  DAG_UNHANDLED_OPTION = moab::MB_UNHANDLED_OPTION,
  DAG_STRUCTURED_MESH = moab::MB_STRUCTURED_MESH,
  DAG_FAILURE = moab::MB_FAILURE
};
#else
//Hard code to take same implicit values as above
enum ErrorCode {
  DAG_SUCCESS = 0,
  DAG_INDEX_OUT_OF_RANGE,
  DAG_TYPE_OUT_OF_RANGE,
  DAG_MEMORY_ALLOCATION_FAILED,
  DAG_ENTITY_NOT_FOUND,
  DAG_MULTIPLE_ENTITIES_FOUND,
  DAG_TAG_NOT_FOUND,
  DAG_FILE_DOES_NOT_EXIST,
  DAG_FILE_WRITE_ERROR,
  DAG_NOT_IMPLEMENTED,
  DAG_ALREADY_ALLOCATED,
  DAG_VARIABLE_DATA_LENGTH,
  DAG_INVALID_SIZE,
  DAG_UNSUPPORTED_OPERATION,
  DAG_UNHANDLED_OPTION,
  DAG_STRUCTURED_MESH,
  DAG_FAILURE
};
#endif

/** Helper class to handle error codes */
class ErrorHandler {
 public:
  ErrorHandler() : _code(DAG_SUCCESS) {};
  ~ErrorHandler() {};
  /** Check the error code for a set error.
   *  Print message if error is detected. */
  virtual void checkSetErr(ErrorCode rval, std::string msg) = 0;

#ifdef MOAB_TYPES_HPP
  /** Overloaded version to cast moab ErrorCode type as DAGMC native ErrorCode */
  void checkSetErr(moab::ErrorCode mbcode, std::string msg) {
    checkSetErr(ErrorCode(mbcode), msg);
  };
#endif
  /** Retrieve the latest error code */
  ErrorCode code() { return _code; };

 protected:
  ErrorCode _code;
};

}

#endif
