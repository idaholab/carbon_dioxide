###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif
# check that MOOSE is available
MOOSE_CONTENT      := $(shell ls $(MOOSE_DIR) 2> /dev/null)
ifeq ($(MOOSE_CONTENT),)
  $(error MOOSE framework does not seem to be available. Make sure that either the submodule is checked out or that your MOOSE_DIR points to the correct location)
endif

FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
###############################################################################
CURRENT_DIR        := $(shell pwd)
CARBON_DIOXIDE_DIR        ?= $(CURRENT_DIR)

# moose submodule status
moose_status := $(shell git -C $(CARBON_DIOXIDE_DIR) submodule status 2>/dev/null | grep moose | cut -c1)
ifneq (,$(findstring +,$(moose_status)))
  ifneq ($(origin CARBON_DIOXIDE_DIR),environment)
    moose_status_msg = "WARNING: Your MOOSE submodule is out of date.\n"
  endif
endif

all:: moose_submodule_status

moose_submodule_status:
	@if [ x$(moose_status_msg) != "x" ]; then echo $(moose_status_msg); fi

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
FLUID_PROPERTIES  := yes
include $(MOOSE_DIR)/modules/modules.mk
###############################################################################


# dep apps
include libSBTL_CarbonDioxide.mk
APPLICATION_DIR    := $(CURRENT_DIR)
APPLICATION_NAME   := carbon_dioxide
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
ADDITIONAL_CPPFLAGS += -Wall -Wextra -Werror
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
#
