EXTERNAL_FILE = ../external/simsopt/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
OUTPUT_DIR = output
VMEC_FILE = $(OUTPUT_DIR)/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc

.PHONY: all comparison clean field.bc booz_xform_field.bc quasi_symmetric.bc test
all: quasi_symmetric.bc

#Alias
quasi_symmetric.bc: $(OUTPUT_DIR)/quasi_symmetric.bc

field.bc: $(OUTPUT_DIR)/field.bc

booz_xform_field.bc: $(OUTPUT_DIR)/booz_xform_field.bc

$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

$(VMEC_FILE): $(EXTERNAL_FILE) $(OUTPUT_DIR)
	cp $< $@

$(OUTPUT_DIR)/quasi_symmetric.bc: $(VMEC_FILE) $(OUTPUT_DIR)
	python make_quasi_symmetric_bc.py $< -1 $@

$(OUTPUT_DIR)/field.bc: $(VMEC_FILE) $(OUTPUT_DIR)
	python run_vmec2boozer.py $< $@

$(OUTPUT_DIR)/booz_xform_field.bc: $(VMEC_FILE) $(OUTPUT_DIR)
	python make_bc_from_booz_xform.py $< $@

comparison: $(OUTPUT_DIR)/field.bc $(OUTPUT_DIR)/booz_xform_field.bc
	python compare_modes.py $(OUTPUT_DIR)

test: booz_xform_field.bc quasi_symmetric.bc $(VMEC_FILE)
	pytest test

clean:
	rm -rf $(OUTPUT_DIR)
