from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd

def get_fits_info(fits_file):
    """從 FITS 檔案提取中心 RA, Dec, 天區半徑及拍攝日期"""
    with fits.open(fits_file) as hdulist:
        header = hdulist[0].header
        wcs = WCS(header)

        # 獲取圖像的中心點和半徑
        img_shape = hdulist[0].data.shape
        center_x, center_y = img_shape[1] // 2, img_shape[0] // 2
        center_coord = wcs.pixel_to_world(center_x, center_y)
        
        corner_coord = wcs.pixel_to_world(0, 0)
        radius = center_coord.separation(corner_coord).deg  # 半徑 (度)

        date_obs = header.get('DATE-OBS', 'N/A')  # 獲取拍攝日期

        return {
            "center_ra": center_coord.ra.deg,
            "center_dec": center_coord.dec.deg,
            "radius_deg": radius,
            "date_obs": date_obs
        }

def read_mpcorb(mpcorb_file):
    """讀取並處理 MPCORB.DAT 檔案的資料，並過濾 Dec 超過範圍的數據"""
    mpc_data = pd.read_csv(mpcorb_file, sep='\s+', skiprows=43, header=None, usecols=[3, 5, 6], 
                           names=["Epoch", "RA_deg", "Dec_deg"], engine='python')
    
    # 過濾掉超出 Dec 範圍的數值
    mpc_data = mpc_data[(mpc_data["Dec_deg"] >= -90) & (mpc_data["Dec_deg"] <= 90)]
    return mpc_data

def find_asteroids_in_fov(mpc_data, center_ra, center_dec, radius_deg):
    """篩選出在給定範圍內的小行星"""
    center_coord = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame='icrs')
    mpc_coords = SkyCoord(ra=mpc_data["RA_deg"].values * u.deg, dec=mpc_data["Dec_deg"].values * u.deg, frame='icrs')
    
    # 計算每個小行星位置與中心點的距離
    sep = center_coord.separation(mpc_coords)
    
    # 篩選出在半徑範圍內的小行星
    within_radius = sep.deg <= radius_deg
    return mpc_data[within_radius]

# 使用範例
fits_file = 'fits/20240811/15h00m00s00d30m00s_060_NoFilt_20240811_000404.fts'
mpcorb_file = 'MPCORB.DAT'

# 提取 FITS 資訊
fits_info = get_fits_info(fits_file)
print(f"拍攝日期: {fits_info['date_obs']}")
print(f"中心 RA: {fits_info['center_ra']} 度, 中心 Dec: {fits_info['center_dec']} 度")
print(f"天區半徑: {fits_info['radius_deg']} 度")

# 讀取 MPCORB 資料
mpc_data = read_mpcorb(mpcorb_file)

# 查找在範圍內的小行星
asteroids_in_fov = find_asteroids_in_fov(mpc_data, fits_info['center_ra'], fits_info['center_dec'], fits_info['radius_deg'])
print(f"在此範圍內的小行星數量: {len(asteroids_in_fov)}")
print(asteroids_in_fov)
