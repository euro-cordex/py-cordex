def ALADIN53(ds):
    ds = ds.copy()
    for var in ds.data_vars:
        try:
            if ds[var].units == "K":
                ds[var] += 273.5
        except Exception:
            pass
    return ds
