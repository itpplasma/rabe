OUTPUT_DIR = output
VMEC_FILE = $(OUTPUT_DIR)/wout_henneberg_qa.nc

.PHONY: all clean quasi_axisymmetric.bc all_modes.bc plot
all: quasi_axisymmetric.bc

#Alias
quasi_axisymmetric.bc: $(OUTPUT_DIR)/quasi_axisymmetric.bc

$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

$(OUTPUT_DIR)/quasi_axisymmetric.bc: $(VMEC_FILE) $(OUTPUT_DIR)
	python ../../python/make_quasi_symmetric_bc.py $< 0 $@

plot: $(OUTPUT_DIR)/quasi_axisymmetric.bc
	python plot_field.py output/quasi_axisymmetric.bc 2
	python plot_field.py output/all_modes.bc 2
	python plot_non_symmetric_modes.py

$(OUTPUT_DIR)/all_modes.bc: $(VMEC_FILE) $(OUTPUT_DIR)
	python ../../python/make_bc_from_booz_xform.py $< $@

clean:
	rm -rf $(OUTPUT_DIR)
