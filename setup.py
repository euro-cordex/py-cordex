#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

# get version
with open("cordex/version.py") as f:
    line = f.readline().strip().replace(" ", "").replace('"', "")
    version = line.split("=")[1]
    __version__ = version


##requirements = [ 'numpy', 'pandas', 'xarray', 'netCDF4' ]
requirements = open("requirements.txt").read().strip().split("\n")

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Lars Buntemeyer",
    author_email='lars.buntemeyer@hzg.de',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Python tools for the Cordex Community.",
    entry_points={
        'console_scripts': [
            'crx-tools=cordex.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='cordex',
    name='py-cordex',
    packages=find_packages(include=['cordex', 'cordex.*', 'cordex.tables.*']),
    package_data={'cordex': ['tables/domains/*.csv', 'tables/data_request/*.csv']},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/euro-cordex/cordex',
    version=__version__,
    zip_safe=False,
)
