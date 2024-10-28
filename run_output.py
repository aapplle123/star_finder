import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

def plot_fits_with_asteroids(fits_file, cat_file):
    """將 FITS 檔案顯示並將 CAT 檔案中小行星的位置標記"""
    # 讀取 FITS 檔案
    with fits.open(fits_file) as hdulist:
        wcs = WCS(hdulist[0].header)
        data = hdulist[0].data
    
    # 顯示 FITS 影像
    plt.figure(figsize=(10, 10))
    plt.imshow(data, cmap='gray', origin='lower', vmin=data.mean() - data.std(), vmax=data.mean() + data.std())
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.title('FITS Image with Asteroids')
    
    # 讀取 CAT 檔案，並確認單位為角度
    asteroids = Table.read(cat_file, format='ascii')
    ra_deg = Angle(asteroids['ALPHA_J2000'], unit=u.deg)
    dec_deg = Angle(asteroids['DELTA_J2000'], unit=u.deg)
    
    # 將 RA/Dec 轉換為 FITS 圖像座標
    asteroid_coords = SkyCoord(ra=ra_deg, dec=dec_deg, frame='icrs')
    x_coords, y_coords = wcs.world_to_pixel(asteroid_coords)
    
    # 在圖上標記小行星位置
    plt.scatter(x_coords, y_coords, edgecolor='red', facecolor='none', s=80, label='Asteroids')
    plt.legend()
    plt.grid(color='white', ls='--')
    
    # 顯示圖像
    plt.show()

# 使用範例
fits_file = 'fits/20240811/15h00m00s00d30m00s_060_NoFilt_20240811_000404.fts'
cat_file = 'asteroids_in_fov.cat'
plot_fits_with_asteroids(fits_file, cat_file)
