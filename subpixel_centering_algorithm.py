#==================================================================================================================
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift, rotate
from skimage import draw

#==================================================================================================================
def subpix_centration_allangles(image, imname, n_angles=36, boxsize=128, **kwargs): 
	'''
		image: image to align
		satradius: radius of saturated data to ignore
		perfect: iterate until the exact center is found to fraction of pix specified perfect=0.1
		boxsize: box size (diameter) within which to measure stddev (default 128 pix diam)
		debug: print debug info messages
	'''

	satradius = kwargs.get( 'satradius', False )
	perfect = kwargs.get( 'perfect', False )
	debug = kwargs.get( 'debug', False )
	header = kwargs.get( 'header', False )

	result, lcx, lcy, oxb, oyb = doRotation(image, n_angles, boxsize, satradius, debug)
	print 'Center of rotational symmetry: ', lcx+oxb, lcy+oyb

	image = np.copy( result )
	if perfect:
		tol = perfect  #tolerance -- size of offset allowed to be ~ zero
		count = 0
		while (abs(oxb) > tol) or (abs(oyb) > tol) :
			result, lcx, lcy, oxb, oyb = doRotation( image, n_angles, boxsize, satradius, debug)
			count += 1
			print count,' = #Additional attempts to get perfect (within tolerance: ', tol ,'pix)'
			print 'Offset from rotational symmetry: ', oxb, oyb

			image = np.copy( result )
		print 'Exiting Perfect while loop, Center of rotational symmetry: ', lcx+oxb, lcy+oyb
	if not header: pf.writeto(imname, result, clobber=True)
	else: pf.writeto(imname, result, header=header, clobber=True)

#==================================================================================================================
def doRotation(image, n_angles, boxsize, satradius, debug):
	targim = np.copy( image )
	target_image = np.copy( image )

	anglearr = ( np.arange( n_angles ) + 0.5 ) * 10.
	oxb_arr = np.zeros( n_angles ) #best offset amount x
	oyb_arr = np.zeros( n_angles ) #best offset amount y

	for k in range(n_angles):
		angle = anglearr[k]
		ref_image = rotate(targim, angle, order=3, reshape=False)

		#subarrays
		box = boxsize #diameter
		r1 = box/2. - box/4.
		r2 = box/2. + box/4.

		ny, nx = target_image.shape
		target = np.copy( target_image[ ny/2. - box/2.: ny/2. + box/2., nx/2. - box/2.: nx/2. + box/2.] )
		refim = np.copy( ref_image[ ny/2. - box/2.: ny/2. + box/2. , nx/2. - box/2.: nx/2. + box/2.] )

		# We carved out a subarray but we'll put it back into the full image at the end
		
		lcx = ( nx - 1 ) / 2. #large array center x
		lcy = ( ny - 1 ) / 2. #large array center y

		# Line up arrays to sub-pixel accuracy
		# subpix.pro
		numy, numx = refim.shape
		nsh = 21 #number of shifts

		tcx = ( numx - 1 ) / 2. #target center x
		tcy = ( numy - 1 ) / 2. #target center y
		rsh = np.zeros (( numy, numx )) #refim shifted
		this_rsh = np.zeros(( numy, numx, nsh, nsh )) #this refim shifted
		this_diff = np.zeros(( numy, numx, nsh, nsh ))
		this_stddev = np.zeros(( nsh, nsh ))

		for i in range( nsh ):  #shift in x direction
			for j in range( nsh ): #shift in y direction
				this_refim = refim
				ox = ( i - 10 ) / 10. #offset x
				oy = ( j - 10 ) / 10. #offset y

				if satradius:
					rr, cc = draw.circle( tcy, tcx, satradius )
					this_refim[ rr,cc ] = -9999999.0
				this_rsh[ :, :, i, j ] = shift( this_refim, (oy, ox), order=3 )
				this_rsh[ :, :, i, j ][ this_rsh[ :, :, i, j ] < -1E6 ] = np.nan
				this_diff[ :, :, i, j ] = this_rsh[ :, :, i, j ] - target
				this_stddev[ i, j ] = np.nanstd( this_diff[ r1:r2, r1:r2, i, j ] )
		bs = np.unravel_index( this_stddev.argmin(), this_stddev.shape )
		oxb = ( ( bs[0] - 10 ) /10. ) / 2. #best offset amount x
		oyb = ( ( bs[1] - 10 ) /10. ) / 2. #best offset amount y
		oxb_arr[k], oyb_arr[k] = oxb, oyb

	if debug:
		plt.plot( anglearr, lcx+oxb_arr, 'k-', lw=1.5 )
		plt.plot( anglearr, lcy+oyb_arr, 'k--', lw=1.5 )
		plt.xlabel( 'Angle of rotational symmetry tested (deg)' )
		plt.ylabel( 'Centroid (pix)' )
		plt.xlim( 0, 360 )
		plt.show()
	oxb = np.mean( oxb_arr ) #Offset x Best
	oyb = np.mean( oyb_arr ) #Offset y Best

	result = shift( target_image, (-oyb, -oxb), order=3 )

	return ( result, lcx, lcy, oxb, oyb )

#==================================================================================================================