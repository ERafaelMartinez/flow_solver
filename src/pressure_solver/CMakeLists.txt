target_sources(${PROJECT_NAME}
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/pressure_solver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/successive_overrelaxation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/red_black_sor.cpp
    ${CMAKE_CURRENT_LIST_DIR}/gauss_seidel.cpp
)
