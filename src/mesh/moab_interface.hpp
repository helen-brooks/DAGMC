#ifndef DAG_MOAB_INTERFACE_HPP
#define DAG_MOAB_INTERFACE_HPP

#include "mesh_interface.hpp"
#include "Moab.hpp"
#include "Error.hpp"

namespace DAGMC {

class ExternalMOAB : public MeshContainer<moab::Interface> {

 public:
  ExternalMOAB(moab::Interface* moabPtrIn) { moabPtr = moabPtrIn; };
  ~ExternalMOAB() {};
  moab::Interface& mesh() override { return *moabPtr; };
  const moab::Interface& const_mesh() const override { return *moabPtr; };
  moab::Interface* mesh_ptr() { return moabPtr; };
  bool isNull() { return (moabPtr == nullptr); };

 private:
  moab::Interface* moabPtr;

};

class ExternalSharedMOAB : public MeshContainer<moab::Interface> {

 public:
  ExternalSharedMOAB(std::shared_ptr<Interface> MBI_shared_ptr) {
    moabPtr = MBI_shared_ptr;
  }
  ~ExternalSharedMOAB() {};

  moab::Interface& mesh()  override { return *moabPtr; };
  const moab::Interface& const_mesh() const override { return *moabPtr; };
  moab::Interface* ptr() override { return moabPtr.get(); };
  std::shared_ptr<moab::Interface> sptr() override { return moabPtr; };
  bool isNull() override { return (moabPtr == nullptr); };

 private:
  std::shared_ptr<moab::Interface> moabPtr;

};


class InternalMOAB : public MeshContainer<moab::Interface> {

 public:
  InternalMOAB() {
    moabPtr = std::make_shared<moab::Core>();
  };
  ~InternalMOAB() {
    moabPtr->delete_mesh();
  };

  moab::Interface& mesh()  override { return *moabPtr; };
  const moab::Interface& const_mesh() const override { return *moabPtr; };
  moab::Interface* ptr() override { return moabPtr.get(); };
  std::shared_ptr<moab::Interface> sptr() override { return moabPtr; };
  bool isNull() override { return (moabPtr == nullptr); };

 private:
  std::shared_ptr<moab::Interface> moabPtr;

};


class MoabInterface : public MeshInterface {

 public:
  MoabInterface(moab::Interface* moabPtrIn);
  MoabInterface(std::shared_ptr<moab::Interface> moabSharedPtrIn);
  ~MoabInterface() {};

  // Read and write to file
  bool load(std::string filename) override;
  bool write(std::string filename) override;

  // Methods for fetching metadata
  bool get_faceting_tol(double& facetingTolerance);
  bool get_group_handles(std::vector<EntityHandle>& group_handles);
  bool get_tag(std::string& tagname, Tag& tag);
  bool get_tag_data(const Tag& tag, const EntityHandle* entityPtr,
                    const int num_handles, void* tag_data);
  bool get_tag_data_arr(const Tag& tag, const EntityHandle* entityPtr,
                        const int num_handles, const void** tag_data, int* len);
  bool get_tag_data_vec(Tag tag, EntityHandle eh, std::vector<std::string>& values);
  bool get_tag_name(Tag tag, EntityHandle eh, std::string& name);
  bool get_group_name(EntityHandle group, std::string& name);
  bool get_entity_sets(EntityHandle group, Range& group_sets);
  bool get_tagged_entity_sets(EntityHandle group, Tag tag, Range& group_sets);
  bool get_tagged_entity_sets(EntityHandle group, std::vector<Tag> tags, Range& group_sets);

  // Methods for setting metadata
  bool set_tag(Tag tag, EntityHandle eh, std::string& new_string);
  bool set_tag_data(Tag tag, const EntityHandle* entityPtr,
                    int num_handles, const void* const tag_data, int len);

  // Retrieve references to moab
  moab::Interface& moab() { return container->mesh(); };
  moab::Interface* moab_ptr() { return container->ptr(); };
  std::shared_ptr<moab::Interface> moab_sptr() { return container->sptr(); };

  // Return error code, but cast as DAGMC error code.
  ErrorCode code() { return ErrorCode(rval); };

  // Get the name tag
  // To-do: what is this for?
  Tag name_tag() { return nameTag; };

 private:

  Tag get_tag(const char* name, int size, TagType store, DataType type,
              const void* def_value = NULL, bool create_if_missing = true);

  // Container for the mesh
  std::shared_ptr<MeshContainer<moab::Interface> > container;

  // Store MOAB return values
  moab::ErrorCode rval;

  Tag nameTag;
  Tag facetingTolTag;

  const int null_delimiter_length;

};

}


#endif
