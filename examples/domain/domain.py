from cordex import domain as dm

# available tables
print(dm.tables())
# available domains names
print(dm.names())
# available domains names
print(dm.names('cordex-core'))
# print cordex core table
print(dm.table('cordex-core'))

# show domains with some dummy data (uses cdo topo)
##for short_name, domain in dm.domains().items():
##    print(domain)
##    domain.to_netcdf(short_name+'.nc', dummy='topo')
