# Overview

NGLView is a miolecular visualization program which runs from a Jupyter notebook environment. This code example is meant to assist using it.

# Required Packages

* mdtraj
    * NGLView can be used with other programs. It just needs a program which can read molecular trajectory/position information. This example uses mdtraj.
* nglview
* ipywidgets=7
* jupyter 

These can be installed using conda:

`conda install -c conda-forge nglview mdtraj ipywidgets=7 jupyter`

# Most Useful Links

* Coloring schemes: https://nglviewer.org/ngl/api/manual/coloring.html 
    * You can also just use common color names i.e. "red", "green", "blue".
* Molecular representations: https://nglviewer.org/ngl/api/manual/molecular-representations.html
* Selection language: https://nglviewer.org/ngl/api/manual/selection-language.html
* Saving a high quality image: https://github.com/nglviewer/nglview/blob/master/examples/notebooks/export_image.ipynb 

# Included Files

* Included in this code example are a PDB file holding a molecular structure, an H5 file holding a molecular dynamics trajectory (see https://mdtraj.org/1.9.4/hdf5_format.html), and a Jupyter Notebook showing how to visualize both.

# Notes
* If you get the error `AttributeError: 'super' object has no attribute '_ipython_display_'
` when you try to import nglview, you probably need to downgrade ipywidgets below version 7.

# Other Useful Links

* NGLView github: https://github.com/nglviewer/nglview
    * README: https://github.com/nglviewer/nglview/blob/master/examples/README.md
    * Thanks to the developers!
* Official documentation: https://nglviewer.org/nglview/latest/
* More examples: https://github.com/nglviewer/nglview/tree/master/examples
