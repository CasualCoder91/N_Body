import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from datetime import datetime #for naming output
import os

import photutils
from photutils import DAOStarFinder, CircularAperture, aperture_photometry
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.background import Background2D, MedianBackground

from astropy.stats import mad_std, sigma_clipped_stats
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch, simple_norm
from astropy.table import Table, QTable
from astropy.nddata import Cutout2D
from astropy import wcs

import config
from database import Database
#from util import magnitude_histogram
from point import Point


def flux_to_mag(flux):
    #http://ircamera.as.arizona.edu/astr_250/Lectures/Lecture_13.htm
    #3880 = 0-magnitude flux in V filter | 640 for K filter
    #Data from http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
    return -2.5 * np.log10(flux/(config.exposure_time*43.6*4000*978e4 ))

def image_segmentation_cut(fits_file):
    n_cuts = 6
    shape = n_pixel/n_cuts # likely n_pixel*2/n_cuts
    w = wcs.WCS(fits_file[0].header)
    stars = QTable()
    for i in np.arange(0,n_cuts,1):
        for j in np.arange(0,n_cuts,1):
            print(i,j)
            position = (i*shape,j*shape)
            cutout = Cutout2D(fits_file[1].data, position, shape-1, wcs=w)
            image = cutout.data[:, :].astype(float)
            stars_segm = image_segmentation(image, False, True)
            if len(stars_segm)>0:
                stars_segm['xcentroid'] = stars_segm['xcentroid']+position[0]-shape/2+2
                stars_segm['ycentroid'] = stars_segm['ycentroid']+position[1]-shape/2+2
                if len(stars)==0:
                    stars = stars_segm
                else:
                    stars = np.append(stars,stars_segm, axis=0)
    stars = stars[[np.isfinite(star['flux']) for star in stars]]
    return stars



def image_segmentation(image, save_color=False,deblend=True):
    #threshold = detect_threshold(image, nsigma=2.)

    box_size = 32 # 8 is too little
    threshold_factor = 50#1.1
    sigma_factor = 5 # FWHM = 3.
    gaussian_size = 8
    contrast = 0.05
    npixel = 9 #für K Filter | J Filter 4

    background2D = Background2D(image,box_size) 
    threshold = threshold_factor * background2D.background_rms 
    data = image - background2D.background  # subtract the background

    sigma = sigma_factor * gaussian_fwhm_to_sigma  
    kernel = Gaussian2DKernel(sigma, x_size=gaussian_size, y_size=gaussian_size)
    kernel.normalize()
    #npixels = minimum amount of connected pixels over threshold to classify as star
    segm = detect_sources(data, threshold, npixels=npixel, filter_kernel=kernel)
    if segm is None:
        return QTable()
    if deblend:
        segm = deblend_sources(image, segm, npixels=npixel,
                                   filter_kernel=kernel, nlevels=32,
                                   contrast=contrast)
    cat = SourceCatalog(data, segm) #from documentation: "data should be background-subtracted for accurate source photometry and properties"
    columns = ['xcentroid', 'ycentroid', 'kron_flux']
    tbl = cat.to_table(columns=columns)
    tbl.rename_column('kron_flux', 'flux')
    #print(tbl)
    #plot
    if(save_img or save_color):
        if(save_color):
            fig = plt.figure(figsize=(n_pixel/960, n_pixel/960), dpi=96)
            cmap = segm.make_cmap(seed=123)
            plt.imshow(segm, origin='upper', cmap=cmap, interpolation='nearest')
            plt.title('Image Segmentation')
            fig.savefig(output_path + "/ISPhotutils_{0}BK_{1}THRE_{2}Kernel{3}{3}_{4}deblended_c{5}.png"
                        .format(box_size,threshold_factor,sigma_factor,gaussian_size,contrast,datetime.now().strftime("%Y_%m_%d_%H%M%S")),dpi=960)
        else:
            positions = np.transpose((tbl['xcentroid'], tbl['ycentroid']))
            apertures = CircularAperture(positions, r=3.)
            fig = plt.figure(figsize=(n_pixel/960, n_pixel/960), dpi=96)
            plt.imshow(image, origin='upper', norm=LogNorm())
            apertures.plot(color='white', lw=0.5, alpha=0.5)
            #fig.savefig(output_path + "/ISPhotutils_{0}BK_{1}THRE_{2}Kernel{3}{3}_{4}deblended.png"
            #            .format(box_size,threshold_factor,sigma_factor,gaussian_size,contrast),dpi=960)
    return tbl

def use_DAOStarFinder(image, save_img=True):

    roundlo = -0.6#0.5
    roundhi = 0.6#0.5
    sharplo = 0.2
    box_size = 32

    mean, median, std = sigma_clipped_stats(image, sigma=3.0)

    #bkg_estimator = MedianBackground()
    #bkg = Background2D(image, box_size, filter_size=(5, 5), bkg_estimator=bkg_estimator)
    #bkg_rms = bkg.background_rms
    threshold = (5. * std) #5.*bkg_sigma

    data = image - median  # subtract the background
    bkg_sigma = mad_std(image)
    daofind = DAOStarFinder(fwhm=3., threshold=threshold)
    mask, center_of_mask = get_masks(image)
    sources = daofind(data, mask=mask)
    if not center_of_mask is None:
        sources = np.hstack([sources, center_of_mask])
    #for col in sources.colnames:
    #    sources[col].info.format = '%.8g'  # for consistent table output
    #print(sources)

    if save_img:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=3.)
        #phot_table = aperture_photometry(image, apertures)
        #for col in phot_table.colnames:
        #    phot_table[col].info.format = '%.8g'  # for consistent table output
        #print(phot_table)
        fig, ax = plt.subplots()
        #fig = plt.figure(figsize=(409.6/96, 409.6/96), dpi=96)
        plt.imshow(image, origin='upper', norm=LogNorm())
        apertures.plot(color='white', lw=0.5, alpha=0.5)
        #for source in sources:
        #    if source["flux"]>200:
        #        width = 0.01*source["flux"]+22
        #        if width > 255:
        #            width = 255
        #        rect = patches.Rectangle((source['xcentroid']-width/2, source['ycentroid']-width/2), width, width, linewidth=1, edgecolor='r', facecolor='none')
        #        ax.add_patch(rect)
        fig.savefig(config.output_path + "/DAOPhotutils_{0}round_{1}sharp_{2}.png".format(roundlo,sharplo,datetime.now().strftime("%Y_%m_%d_%H%M%S")),dpi=960)
    return sources

def get_masks(image):
    mask = np.zeros(image.shape, dtype=bool)

    mean, median, std = sigma_clipped_stats(image, sigma=3.0)

    #bkg_estimator = MedianBackground()
    #bkg = Background2D(image, box_size, filter_size=(5, 5), bkg_estimator=bkg_estimator)
    #bkg_rms = bkg.background_rms
    threshold = (5. * std) #5.*bkg_sigma

    data = image - median  # subtract the background

    bkg_sigma = mad_std(image)
    daofind = DAOStarFinder(fwhm=3., threshold=threshold)
    sources = daofind(data)
    sources.sort(["flux"], reverse=True) #order by flux descending
    center_of_mask = QTable(sources)
    center_of_mask = center_of_mask.remove_rows([0, len(sources)-1])
    #print(len(sources))
    for source in sources:
        if source["flux"]>100:
            width = 0.01*source["flux"]+28
            if width > 255:
                width = 255
            if( mask[int(source['ycentroid']),int(source['xcentroid'])] == False):
                if(center_of_mask is None):
                    center_of_mask = QTable(source)
                else:
                    center_of_mask.add_row(source)
            mask[int(source['ycentroid']-width/2):int(source['ycentroid']+width/2),int(source['xcentroid']-width/2):int(source['xcentroid']+width/2)] = True

    return mask, center_of_mask

#calculates largest distance along x/y between real star and any fake stars.
#returns flux of the real star and found distance
def find_mask_borders(timestep=0):
    hdu = fits.open(output_path +"/scopesim_t"+str(timestep)+".fits")[1]
    image = hdu.data[:, :].astype(float)

    mean, median, std = sigma_clipped_stats(image, sigma=3.0)

    threshold = (5. * std) #5.*bkg_sigma

    data = image - median  # subtract the background

    bkg_sigma = mad_std(image)
    daofind = DAOStarFinder(fwhm=3., threshold=threshold)
    sources = daofind(data)

    max_flux = 0
    max_flux_x = 0
    max_flux_y = 0
    max_flux_dist_x = 0
    max_flux_dist_y = 0

    for source in sources:
        if source["flux"]>max_flux:
            max_flux=source["flux"]
            max_flux_x=source["xcentroid"]
            max_flux_y=source["ycentroid"]

    for source in sources:
        if np.abs(source["xcentroid"]-max_flux_x)>max_flux_dist_x:
            max_flux_dist_x = np.abs(source["xcentroid"]-max_flux_x)
        if np.abs(source["ycentroid"]-max_flux_y)>max_flux_dist_y:
            max_flux_dist_y = np.abs(source["ycentroid"]-max_flux_y)

    if max_flux_dist_x > max_flux_dist_y:
        return max_flux, max_flux_dist_x

    return max_flux, max_flux_dist_y


def pu_all():

    #t0 = np.ndarray((2,),dtype=object)
    #t0[0] = Point(position=np.array([1,2]),velocity = np.array([np.nan,np.nan]),id=1,magnitude=0.5)
    #t0[1] = Point(position=np.array([-4,-5]),velocity = np.array([np.nan,np.nan]),id=2,magnitude=0.1)
    #t1 = np.ndarray((2,),dtype=object)
    #t1[0] = Point(position=np.array([-1.1,-2.2]),velocity = np.array([np.nan,np.nan]),id=-1,magnitude=0.1)
    #t1[1] = Point(position=np.array([1.1,2.2]),velocity = np.array([np.nan,np.nan]),id=-1,magnitude=0.1)
    #generate_velocity_and_index(t0, t1)
    #return

    #print(-2.5* np.log10(17005/(43.6*4000*978e4 )))
    #return
    #db = Database()
    #points_observed = db.select_points(0,True) 
    #points_simulated = db.select_points(0,False) 
    #magnitude_histogram(points_simulated,points_observed)
    #return
    print("Observing stars ...")

    for timestep in [0,1]:
        hdu = fits.open(config.output_path +"/scopesim_t"+str(timestep)+".fits")[1]
        image = hdu.data[:, :].astype(float)

        stars = QTable()
        stars = use_DAOStarFinder(image,False)
        stars = stars[[np.isfinite(star['flux']) for star in stars]]
        # print(stars.info)
        # print("[1] Write found stars to DB!\n[2] no ty")
        print(timestep,": StarFinder done\nWriting stars to DB")

        db = Database()
        origin = config.n_pixel/2. #+0.5 because "For a 2-dimensional array, (x, y) = (0, 0) corresponds to the center of the bottom, leftmost array element. That means the first pixel spans the x and y pixel values from -0.5 to 0.5"
        points = np.ndarray((len(stars),),dtype=object)
        for i, star in enumerate(stars):
            points[i] = Point(position=np.array([config.pixelfactor*(star['xcentroid']-origin),config.pixelfactor*(star['ycentroid']-origin)]),
                                velocity = np.array([np.nan,np.nan]),
                                id=-1,
                                magnitude=flux_to_mag(star['flux']),
                                cluster_id=-1)
        if(timestep>0):
            #points_t0 = db.select_points(timestep-1,True) #get observed stars from previous timestep
            #points_t0, points = generate_velocity_and_index(points_t0,points)
            #db.insert_velocities_HTP(timestep-1,points_t0)
            #db.insert_positions_HTP(timestep,points) #positions with star id = -1 -> asign velocity and id via c++
            db.insert_points(points,timestep)
            #db.update_points(points,timestep)#inserts position and velocity at timestep
        else:
            db.insert_points(points,timestep)

def test():
    timestep=0

    fits_file = fits.open(output_path +"/scopesim_t"+str(timestep)+".fits")
    image = fits_file[1].data[:, :].astype(float)
    #stars = QTable()
    #stars = image_segmentation(image, False, True)
    stars = use_DAOStarFinder(image)
    #stars = image_segmentation_cut(fits_file)

    origin = n_pixel/2. #+0.5 because "For a 2-dimensional array, (x, y) = (0, 0) corresponds to the center of the bottom, leftmost array element. That means the first pixel spans the x and y pixel values from -0.5 to 0.5"
    points = np.ndarray((len(stars),),dtype=object)
    for i, star in enumerate(stars):
        points[i] = Point(position=np.array([pixelfactor*(star['xcentroid']-origin),pixelfactor*(star['ycentroid']-origin)]),
                            velocity = np.array([np.nan,np.nan]),
                            id=-1,
                            magnitude=flux_to_mag(star['flux']),
                            cluster_id=-1)

    db = Database()
    simulated_points = db.select_points(0, False)
    sp_arr = np.vstack(simulated_points[:]).astype(float)
    op_arr = np.vstack(points[:]).astype(float)

    fig = plt.figure()
    fig.set_size_inches(9,9)

    plt.scatter(sp_arr[:,0], sp_arr[:,1], s=1, c='g', marker="o", label='simulated')
    plt.scatter(op_arr[:,0], op_arr[:,1], s=1, c='r', marker="s", label='observed')

    plt.xlabel('ascension [arcsec]', fontsize=16)
    plt.ylabel('declination [arcsec]', fontsize=16)
    plt.legend(loc='upper left')
    plt.show()

if __name__ == '__main__':
    #test()
    pu_all()