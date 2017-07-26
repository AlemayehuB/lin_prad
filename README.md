This is file pertains to the proton radiography reconstruction tool in Python
created by Carlo Graziani, University of Chicago.

Starting files
--------------
Input File: Blob.out
  Note
  ----
    1. Specmesh.py makes specmesh which is the magnetic field and then raytrace.py
     traces all the protons returning a blob.out.
    2. Make an input file from an proton radiography experiment

Progression
-----------
User can either input blob.out like file into either:
  1. reconstruct.py

        Note
        ----
        Takes the input data which contains the initial vector and the final location
        of the proton which the reconstructed magnetic field is calculated. The input
        also contains the Bx and By integral which can be used to get the True Magnetic field

        Returns
        -------
        Fills the second file(blobfield) with the reconstructed magnetic field data and
        true magnetic field data.

        How does it run?
        ---------------
        python reconstruct.py blob.out(experimental data) blobfield.out(Arbitraty name)

  2. image.py

        Note
        ----
        Takes the input data creating various plots using the current integral and
        the final location of the proton. The left panels show the radiograph data,
        the middle panel shows the prediction from the integration, and the right
        panel show the result of subtracting the middle panel from the left panel.

        Returns
        -------
        Outputs the 'Fluence contrast', 'current projection function',
        'Fluence Contrast Difference/Noise','Counts/Bin','Predicted Counts/Bin',
        'Counts Difference/Noise',  'Angular Deviation',and 'Proton Locations' plots
        along with a file that has numerical data that supports the plot.


        How does it run?
        ---------------
        python image.py blob.out(experimental data)

Result
------
  Bplot2.py

    Note
    ----
    Purpose of this script is to compare plots with the True and Reconstructed
    Magnetic field.

    Returns
    -------
    Outputs the projections and deflections of the True and Reconstructed Magnetic Field

    How does it run?
    ---------------
    python Bplot2.py blobfield.out(output of reconstruct.py)
