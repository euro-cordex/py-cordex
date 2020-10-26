"""Example for plotting all orographies with OoPlot by accessing the ESGF database.
"""


import cordex.plot as crxplt

from pyesgf.search import SearchConnection
from pyesgf.logon import LogonManager


# logon to ESGF node
print('logon to ESGF')
lm = LogonManager()
lm.logoff()
lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
print('logged on: {}'.format(lm.is_logged_on()))


def plot_orog(filename, output):
    """plots the orog variable to output file.
    """
    var = 'orog'
    crxplt.contour2(filename, var, output)



# search CORDEX project for REMO2015 fx orog variables
conn = SearchConnection('http://esgf-data.dkrz.de/esg-search', distrib=False)
ctx = conn.new_context(project='CORDEX', experiment='evaluation', time_frequency='fx', rcm_name='REMO2015', variable='orog')
result = ctx.search()

# loop through search results of datasets
for res in result:
    ctx = res.file_context()
    domain = list(ctx.facet_counts['domain'].keys())[0]
    print('domain: {}'.format(domain))
    # the dataset should contains only one files for fx variables
    dataset = ctx.search()
    filename = dataset[0].opendap_url
    print('filename: {}'.format(filename))
    output = '{}.png'.format(domain)
    plot_orog(filename, output)

