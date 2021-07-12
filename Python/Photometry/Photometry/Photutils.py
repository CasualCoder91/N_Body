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

from config import fits_path, output_path, save_img, timestep, pixelfactor,n_pixel
from database import Database

def flux_to_mag(flux):
    #http://ircamera.as.arizona.edu/astr_250/Lectures/Lecture_13.htm
    #3880 = 0-magnitude flux in V filter
    return -2.5 * np.log10(flux/3880)

def image_segmentation(image, save_color=False,deblend=True):
    #threshold = detect_threshold(image, nsigma=2.)

    box_size = 32 # 8 is too little
    threshold_factor = 2
    sigma_factor = 5 # FWHM = 3.
    gaussian_size = 8
    contrast = 0.0005
    npixel = 9 #für K Filter | J Filter 4

    background2D = Background2D(image,box_size) 
    threshold = threshold_factor * background2D.background_rms 
    data = image - background2D.background  # subtract the background

    sigma = sigma_factor * gaussian_fwhm_to_sigma  
    kernel = Gaussian2DKernel(sigma, x_size=gaussian_size, y_size=gaussian_size)
    kernel.normalize()
    #npixels = minimum amount of connected pixels over threshold to classify as star
    segm = detect_sources(data, threshold, npixels=npixel, filter_kernel=kernel)
    if deblend:
        segm = deblend_sources(image, segm, npixels=npixel,
                                   filter_kernel=kernel, nlevels=32,
                                   contrast=contrast)
    cat = SourceCatalog(image, segm)
    columns = ['xcentroid', 'ycentroid', 'kron_flux']
    tbl = cat.to_table(columns=columns)
    #print(tbl)
    #plot
    if(save_img):
        if(save_color):
            fig = plt.figure(figsize=(409.6/96, 409.6/96), dpi=96)
            cmap = segm.make_cmap(seed=123)
            plt.imshow(segm, origin='upper', cmap=cmap, interpolation='nearest')
            plt.title('Image Segmentation')
            fig.savefig(output_path + "/ISPhotutils_{0}BK_{1}THRE_{2}Kernel{3}{3}_{4}deblended_c.png"
                        .format(box_size,threshold_factor,sigma_factor,gaussian_size,contrast),dpi=960)
        else:
            positions = np.transpose((tbl['xcentroid'], tbl['ycentroid']))
            apertures = CircularAperture(positions, r=3.)
            fig = plt.figure(figsize=(409.6/96, 409.6/96), dpi=96)
            plt.imshow(image, origin='upper', norm=LogNorm())
            apertures.plot(color='white', lw=0.5, alpha=0.5)
            fig.savefig(output_path + "/ISPhotutils_{0}BK_{1}THRE_{2}Kernel{3}{3}_{4}deblended.png"
                        .format(box_size,threshold_factor,sigma_factor,gaussian_size,contrast),dpi=960)
    return tbl

#todo: return stars
def use_DAOStarFinder(image):

    mean, median, std = sigma_clipped_stats(image, sigma=3.0)

    bkg_sigma = mad_std(image)
    daofind = DAOStarFinder(fwhm=3., threshold=5.*bkg_sigma)
    sources = daofind(image - median)
    #for col in sources.colnames:
    #    sources[col].info.format = '%.8g'  # for consistent table output
    #print(sources)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    if save_img:
        apertures = CircularAperture(positions, r=3.)
        #phot_table = aperture_photometry(image, apertures)
        #for col in phot_table.colnames:
        #    phot_table[col].info.format = '%.8g'  # for consistent table output
        #print(phot_table)
        fig = plt.figure(figsize=(409.6/96, 409.6/96), dpi=96)
        plt.imshow(image, origin='upper', norm=LogNorm())
        apertures.plot(color='white', lw=0.5, alpha=0.5)
        fig.savefig(output_path + "/DAOPhotutils_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".png",dpi=960)

def main():
    hdu = fits.open(fits_path)[1]
    image = hdu.data[:, :].astype(float)
    stars = image_segmentation(image,save_color=True,deblend=False)
    #use_DAOStarFinder(image)
    #db = Database()
    #origin = n_pixel/2.
    #for star in stars:
    #    db.insert_star(timestep, flux_to_mag(star['kron_flux']), pixelfactor*(star['xcentroid']-origin), pixelfactor*(star['ycentroid']-origin))

if __name__ == '__main__':
    main()