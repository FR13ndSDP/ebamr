AMREX_HOME ?= ../../../../../amrex
TOP := ../

EBASE := NC

DIM := 3

USE_EB := TRUE

LAZY := TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

DEFINES += -DAMREX_CRSEGRNDOMP

include $(TOP)/Source/Make.package
INCLUDE_LOCATIONS += $(TOP)/Source
VPATH_LOCATIONS   += $(TOP)/Source

include $(TOP)/Source/kernel/Make.package
INCLUDE_LOCATIONS += $(TOP)/Source/kernel
VPATH_LOCATIONS   += $(TOP)/Source/kernel

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/EB/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

