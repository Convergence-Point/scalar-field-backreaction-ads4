# Radial scalar field on AdS4. Run `make help` for the targets.
#
#   make              build the animation program and the tests into build/
#   make OPENMP=1     same, with OpenMP (the RK4 loops have parallel pragmas)
#   make check        quick regression test (a few seconds)
#   make test         run the convergence and boundary programs (a few minutes; they print
#                     measured orders, they do not assert anything, see docs/NOTES.md)
#   make animation    run the animation program with the defaults and render animation/wave.gif
#   make convergence  run test_convergence_linear and render convergence/convergence.png
#   make clean        remove build/ (generated data in animation/ and convergence/ is kept)
#
# Layout: src/ is the library, apps/ the animation program, tests/ the measurement programs.
# The programs write to animation/ and convergence/ relative to the current directory,
# so run them from the repository root, e.g. ./build/animation.

CXX      ?= g++
CXXFLAGS ?= -O2
CXXFLAGS += -Wall -Wno-unknown-pragmas -MMD -MP -Isrc
ifeq ($(OPENMP),1)
CXXFLAGS += -fopenmp
endif

BUILD_DIR := build

APP_NAMES  := animation
TEST_NAMES := test_convergence_linear test_convergence_einstein test_convergence_point test_boundary_flat

APPS  := $(addprefix $(BUILD_DIR)/,$(APP_NAMES))
TESTS := $(addprefix $(BUILD_DIR)/tests/,$(TEST_NAMES))

# library objects shared by every program
LIB_OBJS := $(addprefix $(BUILD_DIR)/src/,einstein.o derivative.o evolution.o structure.o)

.PHONY: all clean check test animation convergence help
.DEFAULT_GOAL := all

all: $(APPS) $(TESTS)

# Each program = its own object + the library objects.
$(APPS): $(BUILD_DIR)/%: $(BUILD_DIR)/apps/%.o $(LIB_OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(TESTS): $(BUILD_DIR)/tests/%: $(BUILD_DIR)/tests/%.o $(LIB_OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# ---------------------------------------------------------------------------
# Regression test. Runs in build/check so nothing in animation/ or convergence/
# is overwritten. It checks that
#   1. the linear AdS evolution (mode ads) and its explicitly symmetrized form
#      (mode ads-expanded) agree to 1e-10 in Pi,
#   2. a strong pulse in einstein mode collapses (exit code 3) while a weak one runs to the end.
# ---------------------------------------------------------------------------
check: all
	@rm -rf $(BUILD_DIR)/check && mkdir -p $(BUILD_DIR)/check/animation
	@cd $(BUILD_DIR)/check && ../animation 100 3.1416 20 ads 1 > /dev/null 2>&1 && cp animation/frames.dat ads.dat
	@cd $(BUILD_DIR)/check && ../animation 100 3.1416 20 ads-expanded 1 > /dev/null 2>&1 && cp animation/frames.dat expanded.dat
	@cd $(BUILD_DIR)/check && paste ads.dat expanded.dat | awk 'NF==10 { d = $$2 - $$7; if (d < 0) d = -d; if (d > m) m = d } \
		END { printf "ads vs ads-expanded: max |Pi diff| = %g\n", m; exit (m > 1e-10) }'
	@cd $(BUILD_DIR)/check && ../animation 100 3.1416 20 einstein 1 > /dev/null 2>&1; \
		test $$? -eq 0 && echo "einstein, amplitude 1: runs to the end"
	@cd $(BUILD_DIR)/check && ../animation 100 3.1416 20 einstein 100 > /dev/null 2>&1; \
		test $$? -eq 3 && echo "einstein, amplitude 100: collapse detected (exit 3)"
	@echo "check passed"

# Runs the measurement programs in build/test-run (their output files go there, not into
# animation/ or convergence/). A program that exits with a non-zero code fails the target.
test: all
	@rm -rf $(BUILD_DIR)/test-run && mkdir -p $(BUILD_DIR)/test-run
	@for t in $(TEST_NAMES); do \
		echo "== $$t"; \
		( cd $(BUILD_DIR)/test-run && ../tests/$$t 2>&1 | tail -n 6 ) || exit 1; \
	done
	@echo "all test programs ran"

animation: $(BUILD_DIR)/animation
	./$(BUILD_DIR)/animation
	gnuplot scripts/animation.gp

convergence: $(BUILD_DIR)/tests/test_convergence_linear
	./$(BUILD_DIR)/tests/test_convergence_linear
	gnuplot scripts/convergence.gp

clean:
	rm -rf $(BUILD_DIR)

help:
	@sed -n '1,/^$$/p' Makefile | sed 's/^# \{0,1\}//'

-include $(shell find $(BUILD_DIR) -name '*.d' 2>/dev/null)
