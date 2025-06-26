ifndef SQUID
$(error SQUID is not set)
endif

OUTPUT_DIR = output
VMEC_FILE = wout_squid.nc

all: squid.bc

#Alias
squid.bc: $(OUTPUT_DIR)/squid.bc

$(VMEC_FILE): $(SQUID)
	ln -fs $< $@

$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

plot: squid.bc
	python plot_field.py output/squid.bc 4

$(OUTPUT_DIR)/squid.bc: $(VMEC_FILE) $(OUTPUT_DIR)
	python ../../python/make_bc_from_booz_xform.py $< $@

clean:
	rm -rf $(OUTPUT_DIR)
	rm -rf $(VMEC_FILE)

.PHONY: all clean squid.bc plot
