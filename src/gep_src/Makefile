
BUILDDIR=build

CFLAGS=-Iinclude -Wall -Wextra -pedantic --std=gnu99 -O2 -g
LDFLAGS=-fopenmp

SOURCES=$(shell ls *.c)
SOURCES+=$(shell ls gep_core/*.c)

OBJECTS=$(patsubst %.c, $(BUILDDIR)/%.o, $(notdir $(SOURCES)))

all: $(BUILDDIR) $(BUILDDIR)/libgep.a

$(BUILDDIR)/libgep.a: $(OBJECTS)
	@ar rcs $(BUILDDIR)/libgep.a $(OBJECTS)

$(BUILDDIR)/%.o: %.c
	@gcc $(CFLAGS) $(LDFLAGS) -c $^ -o $@

$(BUILDDIR)/%.o: gep_core/%.c
	@gcc $(CFLAGS) $(LDFLAGS) -c $^ -o $@

$(BUILDDIR):
	@mkdir -p $@

clean:
	rm -rf $(BUILDDIR)
