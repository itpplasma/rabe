VMEC_FILE = ../external/simsopt/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc

.PHONY: all comparison clean
all: field.bc booz_xform_field.bc

field.bc: $(VMEC_FILE)
	python run_vmec2boozer.py $(VMEC_FILE)

booz_xform_field.bc: $(VMEC_FILE)
	python make_bc_from_booz_xform.py $(VMEC_FILE)

quasi_symmetric.bc: $(VMEC_FILE)
	python make_quasi_symmetric_bc.py $(VMEC_FILE) 0

comparison: field.bc booz_xform_field.bc
	python compare_modes.py

clean:
	rm *.bc
