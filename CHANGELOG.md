# GIAS2 Change log

## 0.4.6

This release focuses on the addition of simplex quadratric Lagrange elements and evaluation of arc-lengths to the fieldwork sub-package.

### General
- utility functions for working with simplemeshes
- LineSegment3D setAB method added 
- added version_info variable

### Fieldwork
- new method to create 1d curve elements from sets of nodes
- derivation of quadratric arclength
- test scripts for 1d element arc length evaluation 
- updated hostmeshfit function to working order
- derivation for quadratic simplex elements
- changed topolopy of template 2tri quad patches to be more symmetric
- fixed simplex quadratic lagrange invert mapping 
- completed implementation of simplex_L2_L2 basis

### Image Analysis
- added gaussian_filter method

### Visualisation
- check if scene object is a list of objects when toggling visibility
- draw 1-D mesh using separate curves for each element
- replaced text plotting using mlab.text to mlab.text3d

### OpenSim
- wrapping object class added