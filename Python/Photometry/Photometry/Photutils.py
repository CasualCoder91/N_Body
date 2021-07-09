import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from datetime import datetime #for naming output

import photutils
from photutils import DAOStarFinder, CircularAperture, aperture_photometry
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.background import Background2D

from astropy.stats import mad_std, sigma_clipped_stats
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch, simple_norm

from config import fitsPath, outputPath

def imageSegmentation(image):
    #threshold = detect_threshold(image, nsigma=2.)
    background2D = Background2D(image,32) # 8 is too little
    threshold = 2 * background2D.background_rms #3.0 oder noch mehr?
    data = image - background2D.background  # subtract the background

    sigma = 5 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=8, y_size=8)
    kernel.normalize()
    #npixels = minimum amount of connected pixels over threshold to classify as star
    segm = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)
    segm_deblend = deblend_sources(image, segm, npixels=5,
                               filter_kernel=kernel, nlevels=32,
                               contrast=0.0005)
    cat = SourceCatalog(image, segm_deblend)
    columns = ['xcentroid', 'ycentroid', 'area', 'segment_flux', 'kron_flux']
    tbl = cat.to_table(columns=columns)
    #print(tbl)
    #plot
    if(False):
        fig = plt.figure(figsize=(409.6/96, 409.6/96), dpi=96)
        cmap = segm.make_cmap(seed=123)
        plt.imshow(segm, origin='upper', cmap=cmap, interpolation='nearest')
        plt.title('Image Segmentation')
        #plt.show()
        fig.savefig(outputPath + "/ISPhotutils_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".png",dpi=960)
    if(True):
        positions = np.transpose((tbl['xcentroid'], tbl['ycentroid']))
        apertures = CircularAperture(positions, r=3.)
        fig = plt.figure(figsize=(409.6/96, 409.6/96), dpi=96)
        plt.imshow(image, origin='upper', norm=LogNorm())
        apertures.plot(color='white', lw=0.5, alpha=0.5)
        #for aperture in cat.kron_aperture:
        #    aperture.plot(color='white', lw=0.5)
        fig.savefig(outputPath + "/ISPhotutils_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".png",dpi=960)


def useDAOStarFinder(image):

    mean, median, std = sigma_clipped_stats(image, sigma=3.0)

    bkg_sigma = mad_std(image)
    daofind = DAOStarFinder(fwhm=3., threshold=5.*bkg_sigma)
    sources = daofind(image - median)
    #for col in sources.colnames:
    #    sources[col].info.format = '%.8g'  # for consistent table output
    #print(sources)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=3.)
    #phot_table = aperture_photometry(image, apertures)
    #for col in phot_table.colnames:
    #    phot_table[col].info.format = '%.8g'  # for consistent table output
    #print(phot_table)
    fig = plt.figure(figsize=(409.6/96, 409.6/96), dpi=96)
    plt.imshow(image, origin='upper', norm=LogNorm())
    apertures.plot(color='white', lw=0.5, alpha=0.5)
    fig.savefig(outputPath + "/DAOPhotutils_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".png",dpi=960)

def main():
    hdu = fits.open(fitsPath)[1]
    image = hdu.data[:, :].astype(float)
    imageSegmentation(image)
    #useDAOStarFinder(image)

if __name__ == '__main__':
    main()