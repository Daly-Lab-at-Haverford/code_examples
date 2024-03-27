# Overview

This is a short guide on using matplotlib to make publication quality plots. plot.ipynb is the true "guide", and the other two notebooks have good examples without very much commentary. 

# Guidlines for publication quality plots

* The default size for a plot should be 16 inches wide by 9 inches high. 
    * In some cases, it will make more sense to have a square plot. The most common case is when the x-values and the y-values have some direct comparative relationship.
* Use a large font for axes titles and tick labels. I reccommend 30 for axes titles and 25 for tick labels.
* Make your lines thick (linewidth of 2 or 3) and your markers big (markersize around 15).
* Do NOT use a color for the background of the plot.
* Avoid using a legend unless necessary.
* Do NOT use gridlines.
* **LABEL ALL AXES WITH UNITS!**
* Do NOT connect lines between points unless there is a good scientific or communicative reason. Examples of good reasons:
    * The connecting line(s) represent a real mathematical fit to the data. 
    * The line itself has a physical meaning.
    * The points are very closely spaced, like in a spectrum.
    * The data would be difficult to follow without the line. In this case, use a dotted or dashed line so that it is clear that there is no physical meaning to the line and that it is mearly communicative.
* Avoid giving the plot a title. The plot will appear with a caption in its final location.
* Save your plot as a PDF or in another vector graphics format.
    * If this is somehow not possible, save the figure with a resolution of at least 350 DPI.
* If you're using the same labels, colors, etc in multiple plots, put the labels in a list before 

# Required Packages

* matplotlib
* jupyter 
* notebook
* psi4 is helpful but not required
* scipy is helpful but not required

These can be installed using conda:

`conda install -c conda-forge psi4 numpy scipy matplotlib jupyter notebook`

# Useful Links

* Matplotlib Documentation: https://matplotlib.org/stable/gallery/color/named_colors.html 

# Included Files

* Several Jupyter Notebooks showcasing effective plotting. 

# Possible Improvements

* Include some python scripts for plotting without using Jupyter Notebook.