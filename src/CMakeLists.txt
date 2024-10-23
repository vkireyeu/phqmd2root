if(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

find_package(ZLIB)
include_directories(${CMAKE_SOURCE_DIR}/src)

ROOT_GENERATE_DICTIONARY(G__PHQMDevent PHQMDevent.h LINKDEF LinkDef.h)
add_library(PHQMDevent SHARED PHQMDevent.cxx G__PHQMDevent.cxx)
target_link_libraries(PHQMDevent PUBLIC ${ROOT_LIBRARIES})

string(STRIP ${ROOT_EXE_LINKER_FLAGS} ROOT_EXE_LINKER_FLAGS)

add_executable(phqmd2root phqmd2root.cxx)
target_link_libraries(phqmd2root PUBLIC PHQMDevent ZLIB::ZLIB ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})

add_executable(readphqmd readphqmd.cxx)
target_link_libraries(readphqmd PUBLIC PHQMDevent ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})

add_executable(phqmdtree2ascii phqmdtree2ascii.cxx)
target_link_libraries(phqmdtree2ascii PUBLIC PHQMDevent ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})



install(TARGETS PHQMDevent phqmd2root readphqmd phqmdtree2ascii
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        )

install(FILES
          ${CMAKE_CURRENT_BINARY_DIR}/libPHQMDevent_rdict.pcm
          ${CMAKE_CURRENT_BINARY_DIR}/libPHQMDevent.rootmap
        DESTINATION lib)

install(FILES PHQMDevent.h
        DESTINATION include)

