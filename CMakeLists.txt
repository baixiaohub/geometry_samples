cmake_minimum_required(VERSION 3.18)
project(geometry_samples)

# global configs
set(CMAKE_CXX_STANDARD 14)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# global dirs
set(PROJECT_DEPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/deps)
set(PROJECT_BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(PROJECT_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin/data)

# deps
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)
set(PROJECT_DEPS_LIBS polyscope geometry-central)
set_target_properties(geometry-central PROPERTIES FOLDER "deps")
set_target_properties(polyscope PROPERTIES FOLDER "deps")
set_target_properties(glad PROPERTIES FOLDER "deps")
set_target_properties(imgui PROPERTIES FOLDER "deps")
set_target_properties(stb PROPERTIES FOLDER "deps")
set_target_properties(glfw PROPERTIES FOLDER "deps/glfw")
set_target_properties(update_mappings PROPERTIES FOLDER "deps/glfw")

# targets
add_subdirectory(hello_mesh)
if(MSVC)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT hello_mesh)
endif()
