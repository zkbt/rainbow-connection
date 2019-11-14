'''Tools for resampling array from grid of independent variables to another.'''

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

def binsizes(x):
	'''If x is an array of bin centers, calculate what their sizes are.
		(assumes outermost bins are same size as their neighbors)'''

	binsize = np.zeros_like(x)
	binsize[0:-1] = (x[1:] - x[0:-1])
	binsize[-1] = binsize[-2]
	return binsize

def plotboxy(x, y, **kwargs):
	'''
	Plot with boxes, to show the left and right edges of a box. This is useful,
	for example, to plot flux associated with pixels, in case you are trying to
	do a sub-pixel resample or interpolation or shift.

	(kwargs are passed on to plt.plot)
	'''

	# what are the edges of the bins (making a guess for those on the ends)
	xbinsize = binsizes(x)
	xleft = x - xbinsize/2.0
	xright = x + xbinsize/2.0

	# create a array that doubles up the y values, and interleaves the edges
	plot_x = np.vstack((xleft,xright)).reshape((-1,),order='F')
	plot_y = np.vstack((y,y)).reshape((-1,),order='F')

	# plot those constructed arrays
	plt.plot(plot_x, plot_y, **kwargs)

def fluxconservingresample(xin_unsorted, yin_unsorted, xout, test=False, visualize=False, demo=False,
							treatnanas=0.0):
	'''
	Starting from some initial x and y, resample onto a different grid
	(either higher or lower resolution), while conserving total flux.

	When including the entire range of xin, sum(yout) == sum(yin) should be true.
	'''

	# sort to make sure x is strictly increasing
	s = np.argsort(xin_unsorted)
	xin = xin_unsorted[s]
	yin = yin_unsorted[s]

	# set up the bins, to calculate cumulative distribution of y?
	xinbinsize = binsizes(xin)
	xinleft = xin - xinbinsize/2.0
	xinright = xin + xinbinsize/2.0

	# the first element should be the left edge of the first pixel
	# last element will be right edge of last pixel
	xinforcdf = np.hstack([xinleft, xinright[-1]])

	# to the left of the first pixel, assume flux is zero
	yinforcdf = np.hstack([0, yin])

	# correct for non-finite
	bad = np.isnan(yinforcdf)
	yinforcdf[bad] = treatnanas

	# calculate the cumulative distribution function of the flux (at pixel edge locations)
	cdfin = np.cumsum(yinforcdf)

	# create an interpolator for that
	cdfinterpolator = scipy.interpolate.interp1d(xinforcdf, cdfin,
						kind='linear',
						bounds_error=False,
						fill_value=(0.0, np.sum(yin)))

	# calculate bin edges (of size len(xout)+1)
	xoutbinsize = binsizes(xout)
	xoutleft = xout - xoutbinsize/2.0
	xoutright = xout + xoutbinsize/2.0
	xoutcdf = np.hstack([xoutleft, xoutright[-1]])

	# interpolate the CDF onto those bin edges
	cdfout = cdfinterpolator(xoutcdf)

	# take the derivative of the CDF to get the flux per resampled bin
	# (xout is the center of the bin, and yout is the flux in that bin)
	yout = np.diff(cdfout)

	if visualize:
		fi, (ax_cdf, ax_pdf) = plt.subplots(2,1, sharex=True, figsize=(9,6))
		inkw = dict(color='black', alpha=1, linewidth=3, marker='.', markeredgecolor='none')
		outkw = dict(color='darkorange', alpha=1, linewidth=1, marker='.', markersize=8, markeredgecolor='none')


		legkw = dict(fontsize=10, frameon=False, loc='upper left', bbox_to_anchor=(1, 1))

		# plot the PDFs
		plt.sca(ax_pdf)
		plt.ylabel('Flux per (Original) Pixel')
		plt.xlabel('Pixel')
		# plot the original pixels (in df/dpixel to compare with resampled)
		plotboxy(xin, yin/xinbinsize,
					label='Original Pixels', **inkw)


		# what would a bad interpolation look like?
		badinterpolation = scipy.interpolate.interp1d(xin, yin/xinbinsize,
							kind='linear', bounds_error=False, fill_value=0.0)
		plt.plot(xout, badinterpolation(xout),
						color='cornflowerblue', alpha=1,
						linewidth=1, marker='.',
						markersize=8, markeredgecolor='none',
						label='Silly Simple Interpolation')

		# plot the flux-conserving resampled data (again, in df/d"pixel")
		plt.plot(xout, yout/xoutbinsize, label='Flux-Conserving Interpolation', **outkw)
		plt.legend(**legkw)

		# plot the CDFs
		plt.sca(ax_cdf)
		plt.ylabel('Cumulative Flux (from left)')

		# plot the original CDF
		plt.plot(xinforcdf, cdfin,
					label='Original Pixels', **inkw)

		# plot the interpolated CDF
		plt.plot(xoutcdf, cdfout,
			label='Flux-Conserved Resample', **outkw)
		#plt.legend(**legkw)
		if demo:
			a = raw_input("Pausing a moment to check on interpolation; press return to continue.")

		plt.tight_layout(rect=[0.0, 0.0, 0.67, 1])
		print('{:>6} = {:.5f}'.format('Actual', np.sum(yin)))
		print('{:>6} = {:.5f}'.format('Silly', np.sum(badinterpolation(xout)*xoutbinsize)))
		print('{:>6} = {:.5f}'.format('CDF', np.sum(yout)))


	# return the resampled y-values
	return yout

def bintogrid(x, y, unc=None, newx=None, dx=None, weighting='inversevariance', drop_nans=True):
	'''
	x = the independent variable (wavelength)
	y = the measurement (transit depth)
	unc = the uncertainty on the measurement (sigmas on the transit depths)
	newx = a (linearly) uniformly spaced grid onto which we should bin
	dx = if newx isn't set, then create a grid with this fixed spacing

	Right now, this works only on 1-dimenions arrays,
	but in the long term we should update it to work in
	for arbitrary N-dimensional arrays.

	Arrays must be ordered from decreasing to increasing x.
	'''

	# make up a grid, if one wasn't specified
	if newx is None:
		newx = np.arange(np.min(x), np.max(x) + dx, dx)

	# don't complain about zero-divisions in here (just make uncertainties -> infinity)
	with np.errstate(divide='ignore',invalid='ignore'):

		# resample the sums onto that new grid (in logarithmic space)
		if unc is None:
			weights = np.ones_like(x)
		else:
			if weighting == 'inversevariance':
				weights = 1/unc**2
		numerator = fluxconservingresample(x, y*weights, newx, visualize=False)
		denominator = fluxconservingresample(x, weights, newx, visualize=False)

		# the binned weighted means on the new grid
		newy = numerator/denominator

		# the standard error on the means, for those bins
		newunc = np.sqrt(1/denominator)

	if drop_nans:
		ok = np.isfinite(newy)
	else:
		ok = np.ones_like(newx).astype(np.bool)
	# return binned x, y, unc
	if unc is None:
		return newx[ok], newy[ok]
	else:
		return newx[ok], newy[ok], newunc[ok]


def bintoR(x, y, unc=None, R=50, xlim=None, weighting='inversevariance'):
	'''
	x = the independent variable (wavelength)
	y = the measurement (transit depth)
	unc = the uncertainty on the measurement (sigmas on the transit depths)
	R = the resolution (lambda/dlambda) to bin everything to
	xlim = the limits of the new binned grid
			[if None, will be created from data]

	Right now, this works only on 1-dimenions arrays,
	but in the long term we should update it to work in
	for arbitrary N-dimensional arrays.

	Arrays must be ordered from decreasing to increasing x.
	'''

	# create a new grid of x at the given resolution
	lnx = np.log(x)
	dnewlnx = 1.0/R


	# set the limits of the new xgrid (in log space)
	if xlim is None:
		# use the input grid to set the limits
		lnxbottom, lnxtop = np.min(lnx), np.max(lnx)
	else:
		# use the custom xlim to set the limits
		lnxbottom, lnxtop = xlim

	# create a new, log-uniform grid of x values
	newlnx = np.arange(lnxbottom, lnxtop + dnewlnx, dnewlnx)

	# now do the binning on a uniform grid of lnx

	if unc is None:
		blnx, by = bintogrid(lnx, y, unc, newx=newlnx, weighting=weighting)
		return np.exp(blnx), by
	else:
		blnx, by, bunc = bintogrid(lnx, y, unc, newx=newlnx, weighting=weighting)
		return np.exp(blnx), by, bunc




def testgrid(dx=1.0, N=500, snr=1000, seed=42, irregular=False):
	'''
	a simple test to show how binning to a grid works
	'''


	print('Testing binning {} {} data binned to dx={}\n'.format(N, {True:'irregularly gridded', False:'gridded'}[irregular], dx))


	# set up fake arrays
	np.random.seed(seed)
	if irregular:
		x = np.sort(np.random.uniform(0.5, 25, N))
	else:
		x = np.linspace(0.5, 25, N)
	unc = np.ones_like(x)/snr
	model = np.zeros_like(x) + np.exp(-0.5*(x-5)**2/2**2)*0.01# x*0.0002
	y = np.random.normal(model, unc, N)


	bx, by, bunc = bintogrid(x, y, unc, dx=dx)

	# plot the demonstrations
	plt.figure()
	plt.errorbar(x, y, unc, color='gray', linewidth=0, elinewidth=1, alpha=0.1)
	plt.errorbar(bx, by, bunc, color='black', linewidth=0, elinewidth=3)
	plt.xlabel('x'); plt.ylabel('y')
	plt.title('{} {} data binned to dx={}'.format(N, {True:'irregularly gridded', False:'gridded'}[irregular], dx))
	plt.show()

	# return the results
	return bx, by, bunc


def testR(R=10, N=500, snr=1000, seed=42, irregular=False):
	'''
	a simple test to show how binning to constant R works
	'''

	print('Testing binning {} {} data binned to R={}\n'.format(N, {True:'irregularly gridded', False:'gridded'}[irregular], R))

	# set up fake arrays
	np.random.seed(seed)
	if irregular:
		x = np.sort(np.random.uniform(0.5, 25, N))
	else:
		x = np.linspace(0.5, 25, N)
	unc = np.ones_like(x)/snr
	model = np.zeros_like(x) + np.exp(-0.5*(x-5)**2/2**2)*0.01# x*0.0002
	y = np.random.normal(model, unc, N)

	bx, by, bunc = bintoR(x, y, unc, R=R)

	# plot the demonstrations
	plt.figure()
	plt.errorbar(x, y, unc, color='gray', linewidth=0, elinewidth=1, alpha=0.1)
	plt.errorbar(bx, by, bunc, color='black', linewidth=0, elinewidth=3)
	plt.xscale('log')
	plt.xlabel('x'); plt.ylabel('y')
	plt.title('{} {} data binned to R={}'.format(N, {True:'irregularly gridded', False:'gridded'}[irregular], R))
	plt.show()

	# return the results
	return bx, by, bunc

def testFCR(supersample=True):
	'''this function tests out the resampling code

			supersample=True
				means that there will be multiple new pixels per original pixel

			supersample=False
				means that there will be fewer new pixels than original pixels
	'''

	print("Testing flux-conserving resampling, when {}supersampling the original grid\n".format({True:'', False:'not '}[supersample]))
	xinitial = np.arange(39,47)
	yinitial = np.random.uniform(0.0, 0.1, len(xinitial))
	if supersample:
		xresample = np.linspace(np.min(xinitial) - 1.5, np.max(xinitial) + 1.5, 50)
	else:
		xresample = np.linspace(np.min(xinitial) - 1.5, np.max(xinitial) + 1.5, 5)
	yresample = fluxconservingresample(xinitial, yinitial, xresample, visualize=True)
	plt.show()

def test():
	testFCR()
	testFCR(supersample=False)
	testR(R=10)
	testR(R=10, irregular=True)
	testgrid(dx=0.5)
	testgrid(dx=0.5, irregular=True)
