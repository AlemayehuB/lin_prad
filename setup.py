from setuptools import setup

setup(name='lin-prad',
      version='0.1.0',
      packages=['lin_prad'],
      entry_points={
          'console_scripts': [
              'reconstruct = lin_prad.main:prad_wrap'
          ]
      },
      )
