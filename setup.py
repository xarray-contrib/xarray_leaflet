#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

def get_data_files():
    """Get the data files for the package.
    """
    data_files = [
        ('etc/jupyter/jupyter_server_config.d', ['etc/jupyter/jupyter_server_config.d/xarray_leaflet.json']),
        ('etc/jupyter/jupyter_notebook_config.d', ['etc/jupyter/jupyter_notebook_config.d/xarray_leaflet.json'])
    ]
    return data_files

requirements = [
    'jupyter_server>=0.2.0',
    'rioxarray>=0.0.30',
    'ipyleaflet>=0.13.1',
    'pillow>=7',
    'matplotlib>=3',
    'affine>=2',
    'mercantile>=1'
]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="David Brochart",
    author_email='david.brochart@gmail.com',
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
    description="An xarray extension for map plotting",
    install_requires=requirements,
    license="MIT license",
    long_description="An xarray extension for map plotting",
    include_package_data=True,
    keywords='xarray_leaflet',
    name='xarray_leaflet',
    packages=find_packages(include=['xarray_leaflet', 'xarray_leaflet.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/davidbrochart/xarray_leaflet',
    version='0.1.10',
    zip_safe=False,
    data_files=get_data_files()
)
