import numpy as np
import matplotlib.pyplot as plt
import photutils
from photutils import DAOStarFinder, CircularAperture, aperture_photometry
from photutils.segmentation import detect_threshold, detect_sources
from astropy.stats import mad_std
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch

from config import fitsPath

def imageSegmentation(image):
    threshold = detect_threshold(image, nsigma=2.)

    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(image, threshold, npixels=5, filter_kernel=kernel)

    #plot
    plt.figure(figsize=(12,12))
    cmap = segm.make_cmap(seed=123)
    plt.imshow(segm, origin='lower', cmap=cmap, interpolation='nearest')
    plt.title('Image Segmentation')
    plt.show()


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
    plt.figure(figsize=(12,12))
    plt.imshow(image, cmap='gray_r', origin='upper')
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()

def main():
    hdu = fits.open(fitsPath)[1]
    image = hdu.data[:, :].astype(float)
    imageSegmentation(image)

if __name__ == '__main__':
    main()