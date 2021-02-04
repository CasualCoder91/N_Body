import numpy as np
import photutils
from photutils import DAOStarFinder
from photutils import CircularAperture, aperture_photometry
from astropy.stats import mad_std
from astropy.io import fits
import matplotlib.pyplot as plt

hdu = fits.open(r"E:\Master_Thesis\VS_Project\N_Body\Output\NGC2244\my_image.fits")[1]
image = hdu.data[:, :].astype(float)

bkg_sigma = mad_std(image)
daofind = DAOStarFinder(fwhm=4., threshold=3.*bkg_sigma)
sources = daofind(image)
for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output
print(sources)

positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=4.)
phot_table = aperture_photometry(image, apertures)
for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output
print(phot_table)

plt.imshow(image, cmap='gray_r', origin='upper')
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()
