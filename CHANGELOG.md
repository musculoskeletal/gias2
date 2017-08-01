# GIAS2 Change log

## 0.4.15
- Tri-surface pc registration module and application.
- Colour options for pctraining application.

## 0.4.13
- application scripts

## 0.4.12
- RBF module added the registration subpackage for non-rigid registration.
- Added 3 scripts to the new applications folder
	- rigidreg.py: rigid body registration between point clouds or surface meshes
	- rbfreg.py: RBF-based non-rigid registration between point clouds or surface meshes
	- trainpcashapemodel.py: train a pca-based shape model using a set of correspondent point clouds or surface meshes

## 0.4.11
- updated HJC regression functions to match ISB coordinate system

## 0.4.10
- fixed vector normalisation in fw_segmentation_tools

## 0.4.9
Minor bug fixes.
- tidied up use of common.math.norm in geometric_field
- optional args for array2vtkimage
- fixed spline_tools import

## 0.4.8
- Fixed VTK6 compatibility with image array conversion and binary mask creation

## 0.4.7
- Fixed quadratic simplex element basis functions
- Fixed normalsmoother2 getting edge direction confused with quadratic elements. New method for finding shared edges and their relative directions.

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