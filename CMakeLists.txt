cmake_minimum_required (VERSION 3.8)

project ("Algorithms")

if (NOT TARGET ImaginaCore)
	add_subdirectory ("../Core" Core)
endif ()

add_library (Algorithms MODULE "module.cpp" "Perturbation.cpp" "HarmonicLLA.cpp" "HarmonicMLA.cpp" "MipLA.cpp" "PTWithCompression.cpp")
set_target_properties (Algorithms PROPERTIES PREFIX "")
target_link_libraries (Algorithms ImaginaCore)

set_target_properties (Algorithms PROPERTIES CXX_STANDARD 20)