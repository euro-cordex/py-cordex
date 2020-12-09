

import xesmf as xe


def remap(source, target):    
    regridder = xe.Regridder(source, target, 'bilinear')
    return regridder(source)


