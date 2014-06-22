
BUILDDIR=build

STAT_FILES=$(shell ls ../stats/*)

gen_stats_tables: $(BUILDDIR) scripts/gen_stats_tables.sh $(STAT_FILES)
	@$<

$(BUILDDIR):
	@mkdir -p $@

clean:
	rm -rf $(BUILDDIR)
