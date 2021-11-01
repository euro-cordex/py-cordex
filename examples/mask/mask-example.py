import matplotlib.pyplot as plt

from cordex.regions import _address
from cordex.regions import germany
from cordex import domain as dm


laender_mask = germany.VG2500.regionmask()
laender_geodata = germany.VG2500.geodata()
kreise_mask = germany.VG2500.regionmask("krs")

EUR11 = dm.domain("EUR-11")


fig, ax = plt.subplots(figsize=(13, 10))
fig = laender_mask.plot().get_figure()
fig.savefig("laender.pdf")
fig = kreise_mask.plot().get_figure()
fig.savefig("kreise.pdf")


EUR11_mask = EUR11.gridded_mask(laender_geodata)
EUR11_mask.to_netcdf("EUR-11_laender_mask.nc")


fig, ax = plt.subplots(figsize=(13, 10))
# fig = EUR11_mask.plot(xlim=(-10,0), ylim=(-4,5)).get_figure()
fig = EUR11_mask.plot().get_figure()
fig.savefig("EUR-11_laender_mask.pdf")
