#!/usr/bin/env python

import os

from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(
          name="lin_prad",
          version='1.0.0',

          description='Reconstruction Magnetic Field Tool',
          long_description=open('README.md').read(),

          # Project Homepage
          url='http://flash.uchicago.edu/',

          # Author details
          author='Carlo Graziani and Alemayehu Bogale',
          author_email='carlo@oddjob.uchicago.edu and alemsolobog@uchicago.edu',

          # License
          license='LICENSE.txt',

          # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
          classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science Research',
            'Topic :: Scientific/Engineering :: Physics',

            # Pick your license as you wish (should match "license" above)
             'License :: OSI Approved :: MIT License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 3',
          ],
          # What does it project pertain to?
          keywords='proton radiography',

          ppackages=find_packages(exclude=['contrib', 'docs', 'tests*']),
          install_requires=[
              "numpy >= 1.6",
              "matplotlib >= 2.0",
              "scipy >= 0.19",
            ],
                        # List run-time dependencies here.  These will be installed by pip when
            # your project is installed. For an analysis of "install_requires" vs pip's
            # requirements files see:
            # https://packaging.python.org/en/latest/requirements.html
            #install_requires=['peppercorn'],

            # List additional groups of dependencies here (e.g. development
            # dependencies). You can install these using the following syntax,
            # for example:
            # $ pip install -e .[dev,test]
            #extras_require={
            #    'dev': ['check-manifest'],
            #    'test': ['coverage'],
            #},

            # If there are data files included in your packages that need to be
            # installed, specify them here.  If using Python 2.6 or less, then these
            # have to be included in MANIFEST.in as well.
            #package_data={
            #    'sample': ['package_data.dat'],
            #},

            # Although 'package_data' is the preferred approach, in some case you may
            # need to place data files outside of your packages. See:
            # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
            # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
            #data_files=[('my_data', ['data/data_file'])],

            # To provide executable scripts, use entry points in preference to the
            # "scripts" keyword. Entry points provide cross-platform support and allow
            # pip to create the appropriate form of executable for the target platform.
            entry_points={
                       'console_scripts': ['reconstruct = prad.wrapper:prad'],
                          },
            #},
          )
