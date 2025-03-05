import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

# Load image and catalog
image = fits.getdata('fits/20240811/15h00m00s00d30m00s_060_NoFilt_20240811_000404.fts')
catalog = np.loadtxt('test.cat')

# Plot the image
plt.imshow(image, cmap='gray', origin='lower')

# Overlay detections
plt.scatter(catalog[:, 1], catalog[:, 2], marker='o', color='red', s=10)
plt.show()

# abc