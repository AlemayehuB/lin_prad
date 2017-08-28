#!/usr/bin/env python

import os

from setuptools import setup, find_packages


if __name__ == "__main__":
    setup(name="lin_prad",
          version='1.0.0',
          description='Reconstruction Magnetic Field Tool',
          long_description=open('README').read(),
          license='LICENSE.txt',
          author='Carlo Graziani',
          author_email='carlooddjob.uchicago.edu',
          url='http://flash.uchicago.edu/',
          packages=find_packages(),
          install_requires=[
              "numpy >= 1.6",
            ],
          )
