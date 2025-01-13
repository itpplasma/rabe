VMEC_FILE = ../external/simsopt/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc

.PHONY: all clean
all: field.bc booz_xform_field.bc

field.bc: $(VMEC_FILE)
	python run_vmec2boozer.py $(VMEC_FILE)

booz_xform_field.bc: $(VMEC_FILE)
	python run_booz_xform.py $(VMEC_FILE)

clean:
	rm *.bc