project(TensorVis)

# Minimum version of CMake that can be used
cmake_minimum_required(VERSION 2.6)

# Set up Qt
set(QT_MIN_VERSION "4.5.0")
find_package(Qt4 REQUIRED)

# Set up Avogadro
find_package(Avogadro REQUIRED)
include(${Avogadro_USE_FILE})

# Specify include directories
include_directories(
  # Needed to include Qt moc'd headers
  ${CMAKE_CURRENT_BINARY_DIR}
)

# Create a list of all source files
set( tensorvis_SRCS
     tensorvis.cpp
     #dock.cpp
)

# Create a list of all ui files
set( tensorvis_UIS
     #dock.ui
)

# Create the extension
avogadro_plugin( tensorvis
		 "${tensorvis_SRCS}"
		 "${tensorvis_UIS}")

# Now for the engine:
# Create a list of all source files
set( tensormeshengine_SRCS
     tensormeshengine.cpp
)

avogadro_plugin( tensormeshengine
                 "${tensormeshengine_SRCS}" "")
