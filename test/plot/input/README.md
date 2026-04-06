# field source files

Some original input files for the fields can be found in `test/integration/input`.

Previous input files `.bc` had the wrong sign of the transformation function
`vmns` between toroidal Boozer and geometrical angle for the convension in `neo_field_t`.
The file `landreman_paul_qh_flipped_vmns.bc` has the sign therefore flipped accordingly,
so that it can be compared to the independant calculation of `boozer_field_t`.
