from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd

def get_fits_corners_and_time(fits_file):
    """讀取 FITS 檔案的右上角、左下角座標及拍攝時間"""
    with fits.open(fits_file) as hdulist:
        header = hdulist[0].header
        wcs = WCS(header)
        
        # 擷取拍攝時間
        date_obs = header.get('DATE-OBS', 'N/A')
        
        # 取得圖像的尺寸
        img_shape = hdulist[0].data.shape
        height, width = img_shape
        
        # 定義左下角和右上角的像素座標
        lower_left_pixel = [0, height - 1]
        upper_right_pixel = [width - 1, 0]
        
        # 將像素座標轉換為天球坐標
        lower_left_coord = wcs.pixel_to_world(*lower_left_pixel)
        upper_right_coord = wcs.pixel_to_world(*upper_right_pixel)
        
        return {
            "date_obs": date_obs,
            "lower_left_ra": lower_left_coord.ra.deg,
            "lower_left_dec": lower_left_coord.dec.deg,
            "upper_right_ra": upper_right_coord.ra.deg,
            "upper_right_dec": upper_right_coord.dec.deg
        }

def read_mpcorb(mpcorb_file):
    """讀取並處理 MPCORB.DAT 檔案的資料"""
    mpc_data = pd.read_csv(mpcorb_file, sep='\s+', skiprows=43, header=None, usecols=[3, 5, 6], 
                           names=["Epoch", "RA_deg", "Dec_deg"], engine='python')
    mpc_data = mpc_data[(mpc_data["Dec_deg"] >= -90) & (mpc_data["Dec_deg"] <= 90)]
    return mpc_data

def find_asteroids_in_rectangle(mpc_data, lower_left_ra, lower_left_dec, upper_right_ra, upper_right_dec):
    """篩選出在矩形範圍內的小行星"""
    ra_min, ra_max = sorted([lower_left_ra, upper_right_ra])
    dec_min, dec_max = sorted([lower_left_dec, upper_right_dec])
    in_rectangle = mpc_data[
        (mpc_data["RA_deg"] >= ra_min) & (mpc_data["RA_deg"] <= ra_max) &
        (mpc_data["Dec_deg"] >= dec_min) & (mpc_data["Dec_deg"] <= dec_max)
    ]
    return in_rectangle

def export_to_cat(asteroids, output_cat_file):
    """將小行星資料匯出至 CAT 檔案"""
    with open(output_cat_file, 'w') as file:
        file.write("#   1 NUMBER     Object ID\n")
        file.write("#   2 ALPHA_J2000 Right ascension (J2000) [deg]\n")
        file.write("#   3 DELTA_J2000 Declination (J2000) [deg]\n")
        for i, row in asteroids.iterrows():
            file.write(f"{i+1:<10} {row['RA_deg']:<15.8f} {row['Dec_deg']:<15.8f}\n")

# 使用範例
fits_file = 'fits/20240811/15h00m00s00d30m00s_060_NoFilt_20240811_000404.fts'
mpcorb_file = 'MPCORB.DAT'
output_cat_file = 'asteroids_in_fov.cat'

# 取得 FITS 檔案的邊界和拍攝時間
fits_info = get_fits_corners_and_time(fits_file)
print(f"拍攝日期: {fits_info['date_obs']}")
print(f"左下角座標 - RA: {fits_info['lower_left_ra']} 度, Dec: {fits_info['lower_left_dec']} 度")
print(f"右上角座標 - RA: {fits_info['upper_right_ra']} 度, Dec: {fits_info['upper_right_dec']} 度")

# 讀取 MPCORB 資料
mpc_data = read_mpcorb(mpcorb_file)

# 尋找在 FITS 檔案範圍內的小行星
asteroids_in_rectangle = find_asteroids_in_rectangle(
    mpc_data,
    fits_info['lower_left_ra'],
    fits_info['lower_left_dec'],
    fits_info['upper_right_ra'],
    fits_info['upper_right_dec']
)

# 匯出小行星資料至 CAT 檔案
export_to_cat(asteroids_in_rectangle, output_cat_file)
print(f"在此範圍內的小行星數量: {len(asteroids_in_rectangle)}")
print(f"小行星資料已匯出至 {output_cat_file}")
