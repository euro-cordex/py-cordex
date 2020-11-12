"""this module should help on managing ESGF metadata
"""

import pandas as pd

try:
    from pyesgf.search import SearchConnection
except:
    print('pyesgf client is not installed! please install https://github.com/ESGF/esgf-pyclient to make use of this module.')


DKRZ_HOST = 'esgf-data.dkrz.de'
DKRZ_ESG_URL = 'https://{}/esg-search'.format(DKRZ_HOST)

#cordex  output  EUR-11  GERICS  ECMWF-ERAINT  evaluation  r1i1p1  REMO2015  v1  day  tasmax  v20180813
CORDEX_COLUMNS = ['project_id', 'product', 'domain', 'institute', 'driving_model_id', 'experiment_id', 'member', 'model_id', 'version', 'frequency', 'variable', 'date']
CMIP5_COLUMNS = ['project_id', 'product', 'domain', 'institute', 'driving_model_id', 'experiment_id', 'member', 'model_id', 'version', 'frequency', 'variable', 'date']

COLUMNS  = {'cordex': CORDEX_COLUMNS,
           'cmip5': CMIP5_COLUMNS}


def logon(host=DKRZ_HOST):
    """logon to ESGF Host.
    """
    print('logon to: '.format(host))
    lm = LogonManager()
    lm.logoff()
    lm.logon(hostname=host, interactive=True, bootstrap=True)
    print('logged on: {}'.format(lm.is_logged_on()))
    return lm.is_logged_on()


def connect(url=DKRZ_ESG_URL):
    """establish ESGF search connection.
    """
    return SearchConnection(url, distrib=True)



def search(attrs, conn=None, verbose=False):
    if conn is None:
        return connect().new_context(**attrs).search()
    else:
        return conn.new_context(**attrs).search()

def search_opendap(attrs, conn=None, agg=False, verbose=False):
    """agg=True will return aggregation urls
    instead of list of file urls...
    """
    if conn is None:
        ctx = connect().new_context(**attrs)
    else:
        ctx = conn.new_context(**attrs)
    if verbose: print('Hit Count: {}'.format(ctx.hit_count))
    for facet, item in ctx.facet_counts.items():
        if item: print('{} : {}'.format(facet, item))
    results = ctx.search()
    open_dap = {}
    #print('nr of results: {}'.format(len(results)))
    for result in results:
        sources = []
        if verbose: print(result.dataset_id)
        if agg:
           ctx_results = result.aggregation_context().search()
        else:
           ctx_results = result.file_context().search()
        #print('nr of aggregation results: {}'.format(len(ctx_results)))
        for ctx_result in ctx_results:
            sources.append(ctx_result.opendap_url)
            #print(agg.opendap_url)
        open_dap[result.dataset_id] = sources
    if verbose: print('found {} datasets'.format(len(open_dap.keys())))
    return open_dap



def to_pandas(results):
    """convert a search result to dataframe"""
    data = [res.dataset_id.split('|')[0].split('.') + res.dataset_id.split('|')[1] for res in results]
    #return pd.DataFrame.from_dict(data, orient='index')
    return pd.DataFrame(data, columns = CORDEX_COLUMNS)


