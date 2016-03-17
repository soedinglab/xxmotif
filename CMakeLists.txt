cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(GIMMEmotif)

set(GIMMEMOTIF_MAJOR_VERSION 0)
set(GIMMEMOTIF_MINOR_VERSION 1)
set(GIMMEMOTIF_PATCH_VERSION 1)

set(GIMMEMOTIF_VERSION ${GIMMEMOTIF_MAJOR_VERSION}.${GIMMEMOTIF_MINOR_VERSION}.${GIMMEMOTIF_PATCH_VERSION})

SET (CMAKE_C_COMPILER /usr/bin/cc)
SET (CMKAE_CXX_COMPILER /usr/bin/c++)

set(CMAKE_CXX_FLAGS "-std=c++11 -DLOG_MAX_LEVEL=0 -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -O3 -g3 -pedantic -pedantic-errors -Wall -fmessage-length=0 -fno-strict-aliasing -Wconversion")

set(CMAKE_C_FLAGS "-D__GXX_EXPERIMENTAL_CXX0X__ -O3 -pedantic -pedantic-errors -Wall -fmessage-length=0 -Wconversion")

file(GLOB SOURCES "src/*.cpp" "src/*.h" "src/aminoacids/*.cpp" "src/aminoacids/*.h" "src/elongationPhase/*.cpp" "src/elognationPhase/*.h" "src/em/*.cpp" "src/em/*.h" "src/getopt_pp/*.cpp" "src/getopt_pp/*.h" "src/memoryPool/*.cpp" "src/memoryPool/*.h" "src/nucleotides/*.cpp" "src/nucleotides/*.h" "src/refinementPhase/*.cpp" "src/refinementPhase/*.h" "src/seqFormat/*.cpp" "src/seqFormat/*.h")

add_executable(gimmemotif ${SOURCES})

