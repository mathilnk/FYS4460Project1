TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lattice.cpp \
    atom.cpp \
    facecube.cpp \
    normal.cpp \
    verlet_solver.cpp \
    cellsolver.cpp

HEADERS += \
    lattice.h \
    atom.h \
    facecube.h\
    normal.hpp \
    verlet_solver.h \
    cellsolver.h

