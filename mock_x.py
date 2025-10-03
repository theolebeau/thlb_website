import yt
import pyxsim
import soxs
import sys
import numpy as np
import matplotlib.pyplot as plt

bbox=[[0.4696,0.493,0.4819],[0.4996,0.523,0.5119]]

#fields = [('ramses','x'),('ramses','y'),('ramses','z'),('ramses','Electron_number_density'),('ramses','pressure'),('ramses','temperature')]
#fields = ['x','y','z','Electron_number_density','pressure','temperature']

ds=yt.load("/data/cluster/byopic/SIMS/VirgoClone/output_00251/info_00251.txt", bbox=bbox,default_species_fields="ionized")#, fields=fields)#,fields=["Potential"])
bb=ds.box(left_edge=bbox[0],right_edge=bbox[1])

#print(type(bb))
#print(dir(ds.derived_field_list))
#print(dir(ds.fields.gas))
#sys.exit()

#prj = yt.ProjectionPlot(bb,"z",("ramses", "pressure"),weight_field=("ramses", "cell_volume"),buff_size=(1000, 1000))
#prj.save()

#slc = yt.SlicePlot(bb, "z", ("ramses", "pressure"))#, width=(1.0, "Mpc"))
#slc.show()
#slc.save()

#ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150", default_species_fields="ionized")

#slc = yt.SlicePlot(ds, "z", [("gas", "density"), ("gas", "temperature")], width=(1.0, "Mpc"))
#slc.show()


#source_model = pyxsim.CIESourceModel("spex", 0.05, 11.0, 1000, 0.3, binscale="log")
#xray_fields = source_model.make_source_fields(ds, 0.5, 7.0)
#print(xray_fields)
#sys.exit()
#print(bb["gas", "xray_luminosity_0.5_7.0_keV"])
#print(bb.sum(("gas", "xray_luminosity_0.5_7.0_keV")))

exp_time = (300.0, "ks")  # exposure time
area = (1000.0, "cm**2")  # collecting area
redshift = 0.0038

#n_photons, n_cells = pyxsim.make_photons("sloshing_photons", bb, redshift, area, exp_time, source_model)

#n_events = pyxsim.project_photons("sloshing_photons","sloshing_events","z",(45.0, 30.0),absorb_model="tbabs",nH=0.04)

#soxs.instrument_simulator("sloshing_events.h5","evt.fits",(100.0, "ks"),"chandra_acisi_cy0",[45.0, 30.0],overwrite=True)

#soxs.write_image("evt.fits", "img.fits", emin=0.5, emax=2.0, overwrite=True)

soxs.plot_image("img.fits", stretch="sqrt", cmap="arbre", vmin=0.0, vmax=10.0, width=0.2)
plt.show()

