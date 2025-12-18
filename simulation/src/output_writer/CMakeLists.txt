target_sources(${PROJECT_NAME}
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/output_writer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/output_writer_text.cpp
    ${CMAKE_CURRENT_LIST_DIR}/output_writer_paraview.cpp
    ${CMAKE_CURRENT_LIST_DIR}/output_writer_paraview_parallel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/output_writer_text_parallel.cpp
)
