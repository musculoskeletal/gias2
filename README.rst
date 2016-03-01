GIAS2 (Geometry Image-Analysis Statistics)
==========================================
A Python library for tools used in musculoskeletal modelling. Includes tools for
parametric meshing, registration, image analysis, statistical shape modelling,
and 3-D visualisation using Mayavi.


Dependencies
------------
* scipy
* scikit-learn
* matplotlib


Optional dependencies
---------------------
* VTK and VTK Python bindings (for mesh processing)
* Mayavi (for 3-D visualisation, requires Numpy, VTK, wxPython, configobj)
* PyCSG (for generating constructive solids)
* pydicom (for reading DICOM images)
* Cython (speeds up active shape model and random forest segmentation)


Installation on Linux
---------------------
1. If you would like to use in-built visualisation modules, first install Mayavi for you distribution, else you can skip this step.
    1. Install VTK and VTK python bindings (e.g. through your package manager). VTK 5.10 is the most stable in my experience with Mayavi.
    2. Install mayavi through your package manager (e.g. sudo apt-get install mayavi2) or pip (e.g. pip install --user mayavi)
2. Download the [wheel](https://bitbucket.org/jangle/gias2/downloads/gias2-0.2-py2-none-any.whl) and
    
pip install --user [path/to/wheel]


Installation on Windows
-----------------------
1. The most painless way to install the python dependencies required by GIAS2 is to install the umbrella package [Anaconda](https://www.continuum.io/downloads).
2. If you would like to use in-built visualisation modules, install Mayavi. In you installed Anaconda, from the Anaconda commandline,
        
    conda install mayavi

3. Download the wheel and from the Anaconda commandline
    
    pip install --user [path/to/wheel]


Examples
--------
Example of some the capabilities of GIAS2 can be found in the gias2/examples/ directory. We are working to add more examples.