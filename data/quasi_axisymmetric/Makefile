EXTERNAL_FILE = ../../external/simsopt/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc
OUTPUT_DIR = output
VMEC_FILE = $(OUTPUT_DIR)/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc

.PHONY: all clean quasi_axisymmetric.bc plot
all: quasi_axisymmetric.bc

#Alias
quasi_axisymmetric.bc: $(OUTPUT_DIR)/quasi_axisymmetric.bc

$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

$(VMEC_FILE): $(EXTERNAL_FILE) $(OUTPUT_DIR)
	cp $< $@

$(OUTPUT_DIR)/quasi_axisymmetric.bc: $(VMEC_FILE) $(OUTPUT_DIR)
	python ../../python/make_quasi_symmetric_bc.py $< 0 $@

plot: $(OUTPUT_DIR)/quasi_axisymmetric.bc
	python plot_quasi_axisymmetric.py

clean:
	rm -rf $(OUTPUT_DIR)
