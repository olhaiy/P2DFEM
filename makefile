# Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

-include makefile.local

# Use the MFEM install directory
MFEM_INSTALL_DIR ?= ../mfem-debug
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

EXECUTABLES = batree
EQ_SRC_FILES = SolidConcentration.cpp ElectrolyteConcentration.cpp SolidPotential.cpp ElectrolytePotential.cpp
EQ_INC_FILES = SolidConcentration.hpp ElectrolyteConcentration.hpp SolidPotential.hpp ElectrolytePotential.hpp Equation.hpp
COEFF_INC_FILES = ExchangeCurrentCoefficient.hpp ReactionCurrentCoefficient.hpp OpenCircuitPotentialCoefficient.hpp OverPotentialCoefficient.hpp
SRC_FILES = $(addprefix equations/, $(EQ_SRC_FILES)) P2DOperator.cpp constants.cpp batree.cpp
INC_FILES = $(addprefix equations/, $(EQ_INC_FILES)) $(addprefix coefficients/, $(COEFF_INC_FILES)) parameters.hpp P2DOperator.hpp constants.hpp RegionalCurrent.hpp

.PHONY: all clean

all: $(EXECUTABLES)

# Remove built-in rule
%: %.cpp

batree: $(SRC_FILES) $(INC_FILES) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SRC_FILES) -o $@ $(MFEM_LIBS)

# Replace the default implicit rule for *.cpp files
%: %.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean:
	rm -f *.o *~ $(EXECUTABLES)
	rm -rf *.dSYM *.TVD.*breakpoints
