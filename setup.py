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

          packages=find_packages(exclude=['docs', 'tests*']),
          install_requires=[
              "numpy >= 1.6",
              "matplotlib >= 2.0",
              "scipy >= 0.19",
            ],
        entry_points={
                      'console_scripts': ['reconstruct = prad.wrapper:prad'],
                     },

          )
