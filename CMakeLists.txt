cmake_minimum_required(VERSION 3.16)

project(O2Align LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(O2 REQUIRED)
find_package(OpenMP)

include(AddRootDictionary)

set(O2_ALIGN_PUBLIC_LIBRARIES
    O2::MathUtils
    O2::Steer
    O2::ITSBase
    O2::ITSMFTSimulation
    O2::ITS3Reconstruction
    O2::Framework
    O2::GlobalTrackingWorkflowReaders
    O2::GlobalTrackingWorkflowHelpers
    O2::DataFormatsGlobalTracking
    O2::DetectorsVertexing
    nlohmann_json::nlohmann_json
    GBL::GBL)

add_library(O2Align
            src/DOFSet.cxx
            src/Volume.cxx
            src/Constraint.cxx
            src/Params.cxx
            src/SensorITS.cxx
            src/MisalignmentUtils.cxx # tmp?
            src/AlignmentTypes.cxx
            src/AlignmentSpec.cxx)
#add_library(O2::Align ALIAS O2Align)

set_property(TARGET O2Align PROPERTY EXPORT_NAME O2Align)
set_property(TARGET O2Align PROPERTY OUTPUT_NAME O2Align)

target_include_directories(O2Align
                           PUBLIC
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                             $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                           PRIVATE
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>)

target_link_libraries(O2Align PUBLIC ${O2_ALIGN_PUBLIC_LIBRARIES})

if(OpenMP_CXX_FOUND)
  target_compile_definitions(O2Align PRIVATE WITH_OPENMP)
  target_link_libraries(O2Align PRIVATE OpenMP::OpenMP_CXX)
endif()

add_root_dictionary(O2Align
                    HEADERS include/Align/Params.h
                            include/Align/DOFSet.h
                            include/Align/Volume.h
                            include/Align/Constraint.h
                            include/Align/AlignmentTypes.h
			    LINKDEF src/AlignLinkDef.h)

add_executable(O2exe-alignment-workflow src/alignment-workflow.cxx)
set_property(TARGET O2exe-alignment-workflow
             PROPERTY OUTPUT_NAME o2-its3-alignment-workflow)
target_include_directories(O2exe-alignment-workflow
                           PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(O2exe-alignment-workflow
                      PUBLIC O2Align ROOT::XMLIO)

if(TARGET O2::fpu_support)
  target_link_libraries(O2exe-alignment-workflow PUBLIC O2::fpu_support)
endif()

install(TARGETS O2Align
        EXPORT AlignTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS O2exe-alignment-workflow
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/Align DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/init.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/bin/init.sh
  @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/bin/init.sh
        DESTINATION .
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ
                    WORLD_EXECUTE WORLD_READ)
