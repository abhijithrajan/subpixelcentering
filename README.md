#Algorithm for Subpixel Centering of saturated Point Spread Functions

High contrast imaging searches for exoplanets typically use algorithms that require the knowledge of the star center to subpixel accuracy. For many instruments this problem is complicated by the fact that the star center is either a) saturated, or b) masked by a coronagraph.

The code presented here is converted to Python from an IDL script by Katie Morzinski, described in [Morzinski et al. (2015), ApJ, 815, 108](http://adsabs.harvard.edu/abs/2015ApJ...815..108M)

The algorithm depends on the rotational symmetry of the point spread function. That is to find the center that minimizes the residuals of differenced images over multiple rotation angles. Centers images to subpixel accuracy (for saturated images) by shifting by subpixel (tenth of pixel, from +/- 1 pixel around center). 

The code currently tests rotations over 36 different angles i.e. from 5 to 355 degrees and checks for rotational symmetry in 10 deg. increments and finds the average best centroid -- returns the image centered there.

Limitation of the current code is that improving the centering is done by a brute force grid search and is slow. Future work will focus on improving this.

##Dependencies
numpy
scipy
astropy
python2.7 
scikit-image (if masking the saturated core)
