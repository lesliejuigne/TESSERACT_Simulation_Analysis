add_library(tesssa_processing src/VD_processing.cc)
target_include_directories(tesssa_processing PUBLIC include)
target_link_libraries(tesssa_processing PUBLIC ROOT::Core ROOT::RIO ROOT::Tree)