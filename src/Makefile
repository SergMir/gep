
BUILDDIR=build

CFLAGS=-Igep_src/include -Wall -Wextra -pedantic --std=gnu99 -O2 -g
LDFLAGS=-Wl,--no-as-needed -lm -lGL -lglut -lGLU -fopenmp -Lgep_src/build -lgep

GEP_LIB=gep_src/build/libgep.a
GEP_TEST_BIN=$(BUILDDIR)/gep_test
TEST_FILES_DIR=test_files_gen/build/
COMPARATOR_BIN=$(BUILDDIR)/samples_comparator
RUN_SCRIPT=./scripts/run.sh

usage:
	@echo "Run make with listed arguments:"
	@echo "<> and usage  to show this help"
	@echo "run_single    to run single simple test"
	@echo "run_all_tests to run and collect statistics of all tests"
	@echo "valgrind      to perform memory leaks check with valgrind"
	@echo "clean         to remove all binaries and temporary files"
	@echo "All parameters (such as list of test) could be modified at $(RUN_SCRIPT)"

run_single: build_all
	$(RUN_SCRIPT) origin

run_all_tests: build_all
	$(RUN_SCRIPT) all_tests

valgrind: build_all
	$(RUN_SCRIPT) valgrind

build_all: $(GEP_TEST_BIN) $(COMPARATOR_BIN) gen_test_files

gep_lib:
	@make --directory=gep_src

$(GEP_LIB): gep_lib

$(GEP_TEST_BIN): $(BUILDDIR) $(GEP_LIB) test.c
	@gcc test.c $(CFLAGS) $(LDFLAGS) -o $(GEP_TEST_BIN)

$(COMPARATOR_BIN): samples_comparator.c
	@gcc $^ $(CFLAGS) -lm -o $(COMPARATOR_BIN)

$(BUILDDIR):
	@mkdir -p $@

gen_test_files:
	@make --directory=test_files_gen

clean:
	rm -rf $(BUILDDIR)
	@make --directory=gep_src clean
	@make --directory=test_files_gen clean
