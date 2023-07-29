cmake_minimum_required(VERSION 3.22.1)
project(JSSP)
add_compile_options(-Wno-ignored-attributes)

option(CXX "enable C++ compilation" ON)

if(CXX)
    enable_language(CXX)
endif()
# Set C++ version standard. Using C++20.
set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

# Find GUROBI package.
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

# Set CPLEX installation directory
set(CPLEX_DIR /opt/ibm/ILOG/CPLEX_Studio2211)
# Add CPLEX include directory
include_directories(${CPLEX_DIR}/cplex/include)
include_directories(${CPLEX_DIR}/concert/include)


# Add CPLEX library directory
link_directories(${CPLEX_DIR}/cplex/lib/x86-64_linux/static_pic)
link_directories(${CPLEX_DIR}/concert/lib/x86-64_linux/static_pic)

# Set source files
# Might change to use separate header/source folders.
set (HEADER_FILES
        include/instance.h include/solver.h include/batch_solver.h)

set (SOURCE_FILES
        src/main.cpp src/instance.cpp src/solver.cpp src/batch_solver.cpp)

add_executable(JSSP ${SOURCE_FILES} ${HEADER_FILES})

if(CXX)
    target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
            debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()
# Link CPLEX libraries: MUST BE IN THE RIGHT ORDER!
target_link_libraries(JSSP
        ${CPLEX_DIR}/concert/lib/x86-64_linux/static_pic/libconcert.a
        ilocplex
        ${CPLEX_DIR}/cplex/lib/x86-64_linux/static_pic/libcplex.a
        m
        pthread
        dl
        ${GUROBI_LIBRARY}
        )
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()