"""Create high-res domains from cordex 0.44 domains.
"""


from cordex import domain as dm
import pandas as pd

df = pd.DataFrame()

# create high-res domains from cordex domains
for short_name, domain in dm.domains('cordex').items():
    print('domain: {}'.format(short_name))
    if domain.dlon == 0.44:
        high_res = domain * 0.25
        high_res.short_name = domain.short_name.split('-')[0]+'-11'
        high_res.long_name = domain.long_name
        high_res.region = domain.region
        df = df.append(high_res.to_pandas(), ignore_index=True)


print(dm.table('cordex'))
print(df)
df.to_csv('cordex-high-res.csv', index=False, float_format='%g')
