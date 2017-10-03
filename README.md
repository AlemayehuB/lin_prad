# Linear Proton Radiography Reconstruction

Python package for reconstructing magnetic fields and analyzing proton deflection
data

## Dependencies

* numpy
* matplotlib
* scipy
* pradreader (https://github.com/jtlaune/pradreader)


They can be installed as follows:

```shell
pip install future
pip install numpy matplotlib scipy
pip install git+https://github.com/jtlaune/pradreader
```
For OS X users, it is advised to use Anaconda Python to install `matplotlib` in order to avoid framework errors:

```shell
conda install matplotlib
```


## Installation

This module requires Python 2.7 or 3.5. The latest version can be installed with

```shell
pip install git+https://github.com/AlemayehuB/lin_prad
```

The module can also be installed by

```shell
git clone https://github.com/AlemayehuB/lin_prad
cd lin_prad
python setup.py install
```
## Command Line Tools
### Reconstruct

A command line tool for reconstructing the magnetic field of the data from a proton radiography experiment

Supported input file formats
* Flash
* mitcsv
* carlo

#### Usage
```shell
reconstruct [options] [input file] [file type] [bin length(microns)]
```
##### Options

| Option | Action |
|:-------|--------|
|--tol| The Gauss-Seidel tolerance. DEFAULT:1.0E-04 |
|--iter| The number of Gauss-Seidel iterations. DEFAULT:4000| 

#### Example 
```shell
reconstuct --tol 1.0E-05 --iter 8000 myfile.txt flash4 320
```
This command line script ensures that Gauss-Seidel Tolerance is 1.0E-05 and the number of Gauss-Seidel Iterations 8000 and parses myfile.txt that has a file type of flash4 and a bin length of 320 micron
#### Output

The tool outputs Log Reconstructed Perpendicular Magnetic Field Projection

* If the file has a **carlo** file type then **Path Integrated Magnetci Field Projection**

### Analysis

A command line tool for analysis of a proton radiography experiment

Supported input file formats
* Flash
* mitcsv
* carlo
 
#### Usage
```shell
analysis [input file] [file type]
```
#### Example
```shell
analysis myfile.txt mitcsv
```
This command line parses myfile.txt that has a file type of flash4 and a bin length of 320 micron
#### Output

The tool outputs a Counts/Bin and fluence contrast plot 

* If the file has a **carlo** file type then there is also **Current Projection**, **Predicted Counts/Bin**, and **Noise** for both Counts/Bin and fluence contrast plot.

## Documentation

Full Documentation can be found here
