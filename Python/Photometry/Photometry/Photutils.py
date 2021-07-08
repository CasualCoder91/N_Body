import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from datetime import datetime #for naming output

import photutils
from photutils import DAOStarFinder, CircularAperture, aperture_photometry
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.background import Background2D

from astropy.stats import mad_std
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch, simple_norm

from config import fitsPath, outputPath

def imageSegmentation(image):
    #threshold = detect_threshold(image, nsigma=2.)
    background2D = Background2D(image,32) # 8 is too little
    threshold = 5.0 * background2D.background_rms #3.0 oder noch mehr?
    image -= background2D.background  # subtract the background

    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    #npixels = minimum amount of connected pixels over threshold to classify as star
    segm = detect_sources(image, threshold, npixels=5, filter_kernel=kernel)
    segm_deblend = deblend_sources(image, segm, npixels=5,
                               filter_kernel=kernel, nlevels=32,
                               contrast=0.1)
    cat = SourceCatalog(image, segm_deblend)
    columns = ['xcentroid', 'ycentroid', 'area', 'segment_flux', 'kron_flux']
    tbl = cat.to_table(columns=columns)
    print(tbl)
    #plot
    if(False):
        fig = plt.figure(figsize=(48,48))
        cmap = segm_deblend.make_cmap(seed=123)
        plt.imshow(segm_deblend, origin='upper', cmap=cmap, interpolation='nearest')
        plt.title('Image Segmentation')
        #plt.show()
        fig.savefig(outputPath + "/ISPhotutils_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".png")
    if(False):
        norm = simple_norm(image, 'sqrt')
        fig = plt.figure(figsize=(48,48))
        plt.imshow(image, origin='upper', cmap='Greys_r', norm=LogNorm())
        for aperture in cat.kron_aperture:
            aperture.plot(color='white', lw=1.5)
        fig.savefig(outputPath + "/ISPhotutils_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".png")


def useDAOStarFinder(image):

    image -= np.median(image) #subtract a rough estimate of the background (do not use in combination with detect_threshold) 

    bkg_sigma = mad_std(image)
    daofind = DAOStarFinder(fwhm=4., threshold=3.*bkg_sigma)
    sources = daofind(image)
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
    #print(sources)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.)
    phot_table = aperture_photometry(image, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    #print(phot_table)
    fig = plt.figure(figsize=(12,12))
    plt.imshow(image, cmap='gray_r', origin='upper')
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()

def main():
    hdu = fits.open(fitsPath)[1]
    image = hdu.data[:, :].astype(float)
    imageSegmentation(image)

if __name__ == '__main__':
    main()