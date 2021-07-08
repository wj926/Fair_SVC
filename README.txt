===================================================================
Multi-Basin Support-Based Clustring Toolbox for Matlab
===================================================================

Introduction
------------
The Matlab toolbox MSC-Toolbox implements the support-based clustering methods.
The clustering methods utilize data support functions which user can choose between radial function from the SVDD and variance function of GPR. 

Also user can choose labeling methods, for example 
- complete-graph based labeling method(GP-SC)
- support functions' dynamics based labeling method using stable-equilibrium point(S-MSC)
- transition point based adjacency decision labeling method(T-SMC)
- fast support-based cluster labeling method(F-MSC)
- voronoi cell-based kernel support labeling method(V-MSC) 

After user chooses these types of support function and labeling method, user specifies the options of support function, 
for example, type of kernel function for SVDD or hyper-parameters for GP, and the options for user labeling method. 
Or the program will set methods and options as default values.


Installation
------------
- Create a directory for this toolbox and copy/unzip the toolbox there.
- Go to the toolbox root, ex. (your directory)/msctoolbox. 
- Add path for toolbox directory by running 'msctoolbox.m'. Also you can do this by executing 'msctoolbox' at the command prompt.

Use
---
See 'README.txt' or 'Demo.m' for example.
And see 'msclustering.m' or 'http://slcf.snu.ac.kr/isoftware.html' for explanation for methods and options.
More details are in the each codes of methods.

User can choose or set
- type of support function
- options for support function
- type of labeling method
- options for labeling method
And also can draw plots of the clustering result for 2D data or calculate ARI value.

Version
-------
MSCTOOLBOX version 1.0, 10/10/2014
contact : sujee0524@snu.ac.kr
* The source code is available under the GNU LESSER GENERAL PUBLIC LICENSE, version 2.1.
* This package uses some codes from stprtoolbox.
