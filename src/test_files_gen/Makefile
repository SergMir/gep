
BUILDDIR=build

all: $(BUILDDIR) $(BUILDDIR)/test_sin

$(BUILDDIR)/test_sin: $(BUILDDIR)/test_gen
	@cd $(BUILDDIR) && ./test_gen

$(BUILDDIR)/test_gen: test_gen.c
	@gcc test_gen.c -Wall -pedantic --std=c99 -g -lm -o $(BUILDDIR)/test_gen

$(BUILDDIR):
	@mkdir -p $@

clean:
	rm -rf $(BUILDDIR)