# Copyright 2015 The Edda Authors. All rights reserved.
# Use of this source code is governed by a MIT-style license that can be
# found in the LICENSE file.

# This is only for generating Paraview plugins

cmake_minimum_required (VERSION 3.2)

# List cpp files that will be explosed in plugins
set( PV_SRC_FILES
    eddaUncertainIsocontour.cpp
    eddaRandomSampleField.cpp
    eddaRandomProbeFilter.cpp
    eddaThreshold.cpp
    )

include_directories(${CMAKE_SOURCE_PATH}/src)

if (${EDDA_BUILD_PARAVIEW_PLUGINS})
  ADD_PARAVIEW_PLUGIN(eddaPlugin "1.0"
    SERVER_MANAGER_XML eddaPlugin.xml
    SERVER_MANAGER_SOURCES ${PV_SRC_FILES})
  TARGET_LINK_LIBRARIES(eddaPlugin LINK_PRIVATE edda ${LINK_LIBS})

endif ()


#if(VTK_WRAP_PYTHON)
#  wrap_python(eddaPluginPython "${PV_SRC_FILES}")
#endif()

