# -*- coding: utf-8 -*-
# flake8: noqa
import traceback
import urllib3
import xmltodict



def getxml(url):
    http = urllib3.PoolManager()

    response = http.request('GET', url)
    try:
        data = xmltodict.parse(response.data)
    except:
        print("Failed to parse xml from response (%s)" % traceback.format_exc())
    return data


#def download_cf_table(version='70'):
#    url = 'http://cfconventions.org/Data/cf-standard-names/{}/src/cf-standard-name-table.xml'.format(version)
#    table = requests.get(url, allow_redirects=True)

