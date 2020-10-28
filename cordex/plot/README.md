# Plotting PyRemo and ERA5-Files

This project provides a possibility to create insightful plots with standard python libraries: Xarray and Cartopy.
Used are Data-files from PyRemo with rotated pole coordinates and standard ERA5-files.

## Installation

The packages needed to run the provided scripts are Cartopy, Xarray and NetCDF4. To install them follow the instructions.

```bash
conda install -c conda-forge cartopy
conda install -c conda-forge xarray
pip install netcdf4
pip install ipython
```

Recommended is to create a conda environment. Then follow the instructions to install cordex. And finally, complement the installation with cartopy and ipython.

```bash
conda install -c conda-forge cartopy
pip install ipython
```

###ERA5-Files

The access for the ERA5-Datafiles is described in the GitLab Project.
	https://git.gerics.de/gerics-hiwis/clim-vis/-/tree/master/CDS
Here only one ERA5 file is used to create an animation.


###PyRemo

(This information is taken from the original PyRemo Package. For this project only the Datafiles in OoPlot/examples/data/... are needed. There is no need to proceed with the whole installation to run PyRemo itself.)

You can install the package using pip+git:
    
    pip install git+http://git.gerics.de/python/PyREMO.git

This will install the developmet version by default but you can also specifiy a release:

    pip install git+http://git.gerics.de/python/PyREMO.git@PyRemo-1.2.0
    
Note, that you can define any branch or tag here, but `pip` will not work with
any version lower than 1.2.0.
    
You can also manually clone the package and install it using:

    git clone http://git.gerics.de/python/PyREMO.git
    cd PyREMO
    python setup.py install
    
or install it in developmemt mode:

    pip install -e .
    
which will allow you to edit the source code without reinstalling the package.

# Requirements

To use the full plotting capabilities, you will need to install `pyngl`, please
follwo the instructions here: https://github.com/NCAR/pyngl

I would recommend to first create a conda environment for pyngl and install PyREMO
afterwards in this environment.

##How it works:

The current directory contains following plotting subroutines:

```bash
contour.py
projections.py
stream_vector.py
animation.py
```
Also a subroutine to read the NetCDF provided (read_nc.py) and one to define the plotting characteristics (level_color_label.py).

###Results
The execution routines and outcome of this subroutines in in the following directory:

	examples/plots/

For each plotting subroutine there is an execution routine. The plots are saved in the directory 

	examples/plots/results