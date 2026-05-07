# Generated-style VisIt database plugin build for GR1D restart HDF5 files.
CMAKE_MINIMUM_REQUIRED(VERSION 3.8 FATAL_ERROR)

SET(VISIT_PLUGIN_DIR "$ENV{HOME}/.visit/3.4.2/darwin-arm64/plugins")
INCLUDE("/Applications/VisIt.app/Contents/Resources/3.4.2/darwin-arm64/include/PluginVsInstall.cmake")
INCLUDE("/Applications/VisIt.app/Contents/Resources/3.4.2/darwin-arm64/include/VisItLibraryDependencies.cmake")

PROJECT(GR1DRestartHDF5_database)

SET(COMMON_SOURCES
GR1DRestartHDF5PluginInfo.C
GR1DRestartHDF5CommonPluginInfo.C
)

SET(LIBI_SOURCES
GR1DRestartHDF5PluginInfo.C
)

SET(LIBM_SOURCES
GR1DRestartHDF5MDServerPluginInfo.C
${COMMON_SOURCES}
avtGR1DRestartHDF5FileFormat.C
)

SET(LIBE_SOURCES
GR1DRestartHDF5EnginePluginInfo.C
${COMMON_SOURCES}
avtGR1DRestartHDF5FileFormat.C
)

INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${HDF5_INCLUDE_DIR}
${VISIT_DATABASE_INCLUDES}
)

LINK_DIRECTORIES(${VISIT_LIBRARY_DIR} ${HDF5_LIBRARY_DIR} $ENV{HOME}/.visit/vtk-symlinks)

ADD_LIBRARY(IGR1DRestartHDF5Database ${LIBI_SOURCES})
TARGET_LINK_LIBRARIES(IGR1DRestartHDF5Database visitcommon)
SET(INSTALLTARGETS IGR1DRestartHDF5Database)

IF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)
    ADD_LIBRARY(MGR1DRestartHDF5Database ${LIBM_SOURCES})
    TARGET_LINK_LIBRARIES(MGR1DRestartHDF5Database visitcommon avtdbatts avtdatabase_ser ${HDF5_LIB})
    ADD_TARGET_DEFINITIONS(MGR1DRestartHDF5Database MDSERVER)
    SET(INSTALLTARGETS ${INSTALLTARGETS} MGR1DRestartHDF5Database)
ENDIF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)

ADD_LIBRARY(EGR1DRestartHDF5Database_ser ${LIBE_SOURCES})
TARGET_LINK_LIBRARIES(EGR1DRestartHDF5Database_ser visitcommon avtpipeline_ser avtdatabase_ser ${HDF5_LIB})
ADD_TARGET_DEFINITIONS(EGR1DRestartHDF5Database_ser ENGINE)
SET(INSTALLTARGETS ${INSTALLTARGETS} EGR1DRestartHDF5Database_ser)

IF(VISIT_PARALLEL)
    ADD_PARALLEL_LIBRARY(EGR1DRestartHDF5Database_par ${LIBE_SOURCES})
    TARGET_LINK_LIBRARIES(EGR1DRestartHDF5Database_par visitcommon avtpipeline_par avtdatabase_par ${HDF5_LIB})
    ADD_TARGET_DEFINITIONS(EGR1DRestartHDF5Database_par ENGINE)
    SET(INSTALLTARGETS ${INSTALLTARGETS} EGR1DRestartHDF5Database_par)
ENDIF(VISIT_PARALLEL)

ADD_EXECUTABLE(gr1d_visit_reader_stub gr1d_visit_reader_stub.cpp)
TARGET_INCLUDE_DIRECTORIES(gr1d_visit_reader_stub BEFORE PRIVATE /opt/homebrew/include)
TARGET_LINK_LIBRARIES(gr1d_visit_reader_stub /opt/homebrew/lib/libhdf5.dylib)

VISIT_INSTALL_DATABASE_PLUGINS(${INSTALLTARGETS})
