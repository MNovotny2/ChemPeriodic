FC = mpif90
FFLAGS = -O3 -march=native -ffast-math -funroll-loops -flto -Wall -Wextra
DEBUG = -g -fbacktrace -fbounds-check
PROFILE = -fprofile-generate

# Source file
SRC = ChemPeriodic.f90
EXE = ChemPeriodic

# Check if source exists
ifeq ($(wildcard $(SRC)),)
  $(error Source file $(SRC) not found in $(shell pwd))
endif

.PHONY: all clean debug profile-gen profile-use help

all: $(EXE)

$(EXE): $(SRC)
	@echo "Compiling with $(FC)..."
	$(FC) $(FFLAGS) -o $(EXE) $(SRC)
	@echo "Build successful: ./$(EXE)"

debug: $(SRC)
	@echo "Compiling debug version..."
	$(FC) $(DEBUG) -O1 -o $(EXE)_debug $(SRC)
	@echo "Debug build successful: ./$(EXE)_debug"

profile-gen: $(SRC)
	@echo "Compiling for profiling (generation)..."
	$(FC) $(FFLAGS) $(PROFILE) -o $(EXE)_profile $(SRC)
	@echo "Profile generation build successful. Run ./$(EXE)_profile with test data."

profile-use: $(SRC)
	@echo "Compiling with profile data..."
	$(FC) $(FFLAGS) -fprofile-use -fprofile-correction -o $(EXE) $(SRC)
	@echo "Profile-optimized build successful: ./$(EXE)"

clean:
	@echo "Cleaning build artifacts..."
	rm -f $(EXE) $(EXE)_debug $(EXE)_profile *.o *.mod *.gcda *.gcno
	@echo "Clean complete"

help:
	@echo "Available targets:"
	@echo "  make              - Build optimized binary"
	@echo "  make debug        - Build with debug symbols (-g)"
	@echo "  make profile-gen  - Build for PGO collection"
	@echo "  make profile-use  - Rebuild using PGO data"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make help         - Show this message"
