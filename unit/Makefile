###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# HERD_TRUNK_DIR   - Location of the HERD repository
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
CARBON_DIOXIDE_DIR			?= $(shell dirname `pwd`)
MOOSE_SUBMODULE    := $(CARBON_DIOXIDE_DIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname $(CARBON_DIOXIDE_DIR))/moose
endif
# check that MOOSE is available
MOOSE_CONTENT      := $(shell ls $(MOOSE_DIR) 2> /dev/null)
ifeq ($(MOOSE_CONTENT),)
  $(error MOOSE framework does not seem to be available. Make sure that either the submodule is checked out or that your MOOSE_DIR points to the correct location)
endif

MOOSE_DIR          ?= $(CARBON_DIOXIDE_DIR)/moose
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
###############################################################################
CURRENT_DIR        := $(shell pwd)

# Extra stuff for GTEST
ADDITIONAL_INCLUDES += -I$(FRAMEWORK_DIR)/contrib/gtest
ADDITIONAL_LIBS     += $(FRAMEWORK_DIR)/contrib/gtest/libgtest.la

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk
################################## MODULES ####################################
FLUID_PROPERTIES  := yes
include           $(MOOSE_DIR)/modules/modules.mk
###############################################################################

include $(CARBON_DIOXIDE_DIR)/libSBTL_CarbonDioxide.mk

# dep apps
APPLICATION_DIR    := $(CARBON_DIOXIDE_DIR)
APPLICATION_NAME   := carbon_dioxide
include            $(FRAMEWORK_DIR)/app.mk

APPLICATION_DIR    := $(CURRENT_DIR)
APPLICATION_NAME   := carbon_dioxide-unit
BUILD_EXEC         := yes
app_BASE_DIR       :=      # Intentionally blank
DEP_APPS           ?= $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include $(FRAMEWORK_DIR)/app.mk

# Find all the carbon dioxide unit test source files and include their dependencies.
carbon_dioxide-unit_srcfiles := $(shell find $(CURRENT_DIR) -name "*.C")
carbon_dioxide-unit_deps := $(patsubst %.C, %.$(obj-suffix).d, $(carbon_dioxide-unit_srcfiles))
-include $(carbon_dioxide-unit_deps)

###############################################################################
# Additional special case targets should be added here
