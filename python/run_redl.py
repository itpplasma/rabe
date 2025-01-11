import os
from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer

file_name = os.path.join("external/simsopt/tests/test_files/",
                         "wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc")
vmec = Vmec(file_name)
boozer = Boozer(vmec)
