from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np

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
            "date_obs": date_obs,
            "wcs": wcs
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

def export_to_cat(fits_info, asteroids, output_file):
    """將找到的小行星匯出為 .cat 檔案，包含特定的欄位結構"""
    wcs = fits_info["wcs"]
    with open(output_file, 'w') as f:
        # 編寫標題行
        f.write("#   1 NUMBER                 Running object number\n")
        f.write("#   2 X_IMAGE                Object position along x                                    [pixel]\n")
        f.write("#   3 Y_IMAGE                Object position along y                                    [pixel]\n")
        f.write("#   4 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]\n")
        f.write("#   5 DELTA_J2000            Declination of barycenter (J2000)                          [deg]\n")
        f.write("#   6 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]\n")
        f.write("#   7 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]\n")
        f.write("#   8 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]\n")
        f.write("#   9 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]\n")
        f.write("#  10 FLAGS                  Extraction flags\n")
        
        # 生成每個小行星的數據行
        for idx, row in asteroids.iterrows():
            # 轉換 RA, Dec 為圖像內的 pixel 座標
            sky_coord = SkyCoord(ra=row['RA_deg'] * u.deg, dec=row['Dec_deg'] * u.deg, frame='icrs')
            x_img, y_img = wcs.world_to_pixel(sky_coord)
            
            # 隨機生成其他所需數據
            mag_auto = np.random.uniform(-12, -8)   # 隨機亮度值
            magerr_auto = np.random.uniform(0.01, 0.5)  # 隨機亮度誤差
            flux_auto = np.random.uniform(1000, 10000)  # 隨機光通量
            fluxerr_auto = np.random.uniform(100, 1200)  # 隨機光通量誤差
            flags = np.random.choice([0, 16, 24])  # 隨機標誌值
            
            # 寫入 .cat 格式行
            f.write(f"{idx + 1:>9} {x_img:>10.4f} {y_img:>10.4f} {row['RA_deg']:>12.7f} {row['Dec_deg']:>12.7f} "
                    f"{mag_auto:>8.4f} {magerr_auto:>8.4f} {flux_auto:>12.2f} {fluxerr_auto:>12.2f} {flags:>3}\n")

# 使用範例
fits_file = 'fits/20240811/15h00m00s00d30m00s_060_NoFilt_20240811_000404.fts'
mpcorb_file = 'MPCORB.DAT'
output_cat_file = 'output_asteroids.cat'  # 匯出文件名稱

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

# 將結果匯出為 .cat 檔案
export_to_cat(fits_info, asteroids_in_fov, output_cat_file)
print(f"已匯出小行星資料到 '{output_cat_file}'")
