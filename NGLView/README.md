# Overview

NGLView is a miolecular visualization program which runs from a Jupyter notebook environment. This code example is meant to assist using it.

# Required Packages

* mdtraj
    * NGLView can be used with other programs. It just needs a program which can read molecular trajectory/position information. This example uses mdtraj.
* nglview
* jupyter 
* notebook

These can be installed using conda:

`conda install -c conda-forge nglview mdtraj jupyter notebook`

# Most Useful Links

* Coloring schemes: https://nglviewer.org/ngl/api/manual/coloring.html 
    * You can also just use common color names i.e. "red", "green", "blue".
* Molecular representations: https://nglviewer.org/ngl/api/manual/molecular-representations.html
* Selection language: https://nglviewer.org/ngl/api/manual/selection-language.html
* Saving a high quality image: https://github.com/nglviewer/nglview/blob/master/examples/notebooks/export_image.ipynb 

# Included Files

* Included in this code example are a PDB file for dinitrobenzine.pdb, a PDB file for a solvated acyl carrier protein structure, and a Jupyter Notebook showing how to visualize both.

# Notes
* If you get the error `AttributeError: 'super' object has no attribute '_ipython_display_'
` when you try to import nglview, you probably should reinstall nglview, jupyter, and notebook

# Other Useful Links

* NGLView github: https://github.com/nglviewer/nglview
    * README: https://github.com/nglviewer/nglview/blob/master/examples/README.md
    * Thanks to the developers!
* Official documentation: https://nglviewer.org/nglview/latest/
* More examples: https://github.com/nglviewer/nglview/tree/master/examples

# Citation

Hai Nguyen, David A Case, Alexander S Rose, NGLview–interactive molecular graphics for Jupyter notebooks, Bioinformatics, Volume 34, Issue 7, April 2018, Pages 1241–1242, https://doi.org/10.1093/bioinformatics/btx789