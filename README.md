> It's a very good start to documentation for release! I have a few comments and we'll iterate. -Scott

# Linear Proton Radiography Reconstruction

> A longer description here would be helpful! One or two paragraphs.
> Describe what the use of this package is; e.g., taking a radiograph and reconstructing magnetic field map.
> Cite Carlo's paper as the basis for this code: https://arxiv.org/abs/1603.08617

Python package for reconstructing magnetic fields and analyzing proton deflection
data based on proton radiographyflux image

## Dependencies
This module requires **Python 2.7** or **3.5**. Installation requires **git**.

**OS X users:** Prior to installing dependencies, ensure an adequate Python installation (non-Apple) by following [this guide](https://matplotlib.org/faq/installing_faq.html#osx-notes).

The following Python packages are required:
* numpy (Scientific computing)
* matplotlib (Plotting)
* scipy (Scientific computing)
* future (Cross-compatibility between Python2 and Python3)
* pradreader (https://github.com/jtlaune/pradreader) (Reading various proton radiograph file formats)

On most systems (see above note for OS X), they can be installed using Python's [PIP package manager](https://packaging.python.org/tutorials/installing-packages/) as follows:

```shell
pip install future
pip install numpy matplotlib scipy
pip install git+https://github.com/jtlaune/pradreader
```

Depending on how Python was installed on your system, `pip` may require Administrative or `sudo` privileges.

## Installation
After installing dependencies, install the latest version of **lin_prad** by:

```shell
pip install git+https://github.com/AlemayehuB/lin_prad
```

The module can also be installed by:

```shell
git clone https://github.com/AlemayehuB/lin_prad
cd lin_prad
python setup.py install
```
## Requirements
> This description should clarify that pradreader makes the file below.

An intermediate file that contains the variables such as:
* Distance from proton source to the interaction region(cm), s2r_cm
* Distance from proton source to the screen(cm), s2d_cm 
* Proton Kinetic Energies (MeV), Ep_MeV
* flux image which is a matrix with number of protons per bin of the screen which is dependent on the inputted bin length, flux
* flux reference which is the flux image if there were no interaction region


## Command Line Tools
### Reconstruct

A command line tool for reconstructing the magnetic field of the data from a proton radiography experiment

> Does THIS code support these file formats, or does pradreader? Maybe we could to say something like: "Supported file formats include any that pradreader supports, including radiographs generated from FLASH simulations and MIT's CR39 proton radiography analysis"

Supported input file formats
* Flash
* mitcsv
* carlo

> Outsiders don't know what a "carlo" file type is, nor a "mitcsv"; these names were made up for internal use in pradreader

#### Usage
```shell
reconstruct [options] [input file] [file type] [bin length(microns)]
```
##### Options

| Option | Action |
|:-------|--------|
|--tol| The Gauss-Seidel tolerance. DEFAULT:1.0E-04 |
|--iter| The number of Gauss-Seidel iterations. DEFAULT:4000|


> What do these do? Why would I change these values? What happens if I don't set them? Add a qualitative description of how to use these values (e.g. "Modifying the Gauss-Seidel tolerance changes the threshold for convergence; increasing may speed up convergence but result in lower quality", and cite the section of Carlo's paper that tells them more.

#### Example 
```shell
reconstuct --tol 1.0E-05 --iter 8000 myfile.txt flash4 320
```
This command line script ensures that Gauss-Seidel Tolerance is 1.0E-05 and the number of Gauss-Seidel Iterations 8000 and parses myfile.txt that has a file type of flash4 and a bin length of 320 micron
#### Output

The tool outputs Log Reconstructed Perpendicular Magnetic Field Projection

* If the file has a **carlo** file type then **Path Integrated Magnetic Field Projection**
> What does this mean?

### Analysis

A command line tool for analysis of a proton radiography experiment
> Again, more description here. Why would an outsider want to use this?

Supported input file formats
* Flash
* mitcsv
* carlo
> Same comments as in the other file format section above

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


