"""this module should help on managing ESGF metadata
"""

import pandas as pd

try:
    from tqdm import tqdm
except Exception:

    def tqdm(x):
        return x


try:
    from pyesgf.search import SearchConnection
except Exception:
    print(
        "pyesgf client is not installed! please install https://github.com/ESGF/esgf-pyclient to make use of this module."
    )


DKRZ_HOST = "esgf-data.dkrz.de"
DKRZ_ESG_URL = "https://{}/esg-search".format(DKRZ_HOST)

# cordex  output  EUR-11  GERICS  ECMWF-ERAINT  evaluation  r1i1p1  REMO2015  v1  day  tasmax  v20180813
CORDEX_COLUMNS = [
    "project_id",
    "product",
    "domain",
    "institute",
    "driving_model_id",
    "experiment_id",
    "member",
    "model_id",
    "version",
    "frequency",
    "variable",
    "date",
    "host",
]
CMIP5_COLUMNS = [
    "project_id",
    "product",
    "domain",
    "institute",
    "driving_model_id",
    "experiment_id",
    "member",
    "model_id",
    "version",
    "frequency",
    "variable",
    "date",
]
CMIP6_COLUMNS = [
    "project_id",
    "product",
    "domain",
    "institute",
    "driving_model_id",
    "experiment_id",
    "member",
    "model_id",
    "version",
    "frequency",
    "variable",
    "date",
]


def logon(host=DKRZ_HOST):
    """logon to ESGF Host."""
    from pyesgf.logon import LogonManager

    print("logon to: {}".format(host))
    lm = LogonManager()
    lm.logoff()
    lm.logon(hostname=host, interactive=True, bootstrap=True)
    print("logged on: {}".format(lm.is_logged_on()))
    return lm.is_logged_on()


def connect(url=DKRZ_ESG_URL):
    """establish ESGF search connection."""
    return SearchConnection(url, distrib=True)


def context(attrs, conn=None, verbose=False):
    """get an ESGF search context."""
    if conn is None:
        ctx = connect().new_context(**attrs)
    else:
        ctx = conn.new_context(**attrs)
    if verbose:
        print("Hit Count: {}".format(ctx.hit_count))
    return ctx


def datasets(attrs, conn=None, verbose=False):
    """returns dataset results"""
    return context(attrs, conn, verbose).search()


def dataset_dict(attrs, conn=None, verbose=False):
    return {ds.dataset_id: ds for ds in tqdm(datasets(attrs, conn, verbose))}


def get_opendap_urls(attrs, conn=None, agg=False, verbose=False):
    """agg=True will return aggregation urls
    instead of list of file urls...
    """
    ctx = context(attrs, conn, verbose)
    # for facet, item in ctx.facet_counts.items():
    #    if item: print('{} : {}'.format(facet, item))
    results = ctx.search()
    open_dap = {}
    for result in tqdm(results):
        sources = []
        if verbose:
            print(result.dataset_id)
        if agg:
            ctx_results = result.aggregation_context().search()
        else:
            ctx_results = result.file_context().search()
        for ctx_result in ctx_results:
            sources.append(ctx_result.opendap_url)
        open_dap[result.dataset_id] = sources
    if verbose:
        print("found {} datasets".format(len(open_dap.keys())))
    return open_dap


def _split_dataset_id(dataset_id):
    """splits a dataset id into list of values."""
    return dataset_id.split("|")[0].split(".") + [(dataset_id.split("|")[1])]


def _get_columns(columns=None):
    if columns == "CORDEX":
        return CORDEX_COLUMNS
    elif columns == "CMIP5":
        return CMIP5_COLUMNS
    elif columns == "CMIP6":
        return CMIP6_COLUMNS
    else:
        return columns


def to_pandas(results, columns=None):
    """convert a search result to dataframe"""
    data = [_split_dataset_id(res.dataset_id) for res in tqdm(results)]
    return pd.DataFrame(data, columns=_get_columns(columns))
