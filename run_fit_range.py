from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def analyze_fits(fits_file):
    # 開啟 FITS 檔案並提取頭資訊
    with fits.open(fits_file) as hdulist:
        header = hdulist[0].header
        wcs = WCS(header)

        # 獲取 FITS 拍攝日期和時間（假設關鍵字為 'DATE-OBS'）
        date_obs = header.get('DATE-OBS', 'N/A')

        # 獲取圖像的尺寸
        img_shape = hdulist[0].data.shape
        center_x, center_y = img_shape[1] // 2, img_shape[0] // 2

        # 轉換圖像中心的像素坐標為天球座標
        center_coord = wcs.pixel_to_world(center_x, center_y)

        # 計算圖像對角像素距離，估計天區半徑
        corner_coord = wcs.pixel_to_world(0, 0)
        radius = center_coord.separation(corner_coord).deg  # 單位：度

        return {
            "date_obs": date_obs,
            "center_ra": center_coord.ra.deg,
            "center_dec": center_coord.dec.deg,
            "radius_deg": radius
        }

# 使用範例
fits_file = 'fits/20240811/15h00m00s00d30m00s_060_NoFilt_20240811_000404.fts'
result = analyze_fits(fits_file)
print(f"拍攝日期與時間: {result['date_obs']}")
print(f"中心 RA: {result['center_ra']} 度")
print(f"中心 Dec: {result['center_dec']} 度")
print(f"天區半徑: {result['radius_deg']} 度")
