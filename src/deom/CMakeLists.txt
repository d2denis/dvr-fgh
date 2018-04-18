file (GLOB cur *.cpp)
add_library (fgh ${cur})
add_dependencies(fgh json11)
