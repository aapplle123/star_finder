from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd

# 讀取 MPCORB.DAT 資料
def load_mpc_data(file_path):
    # 使用精確分隔符來解析 MPCORB 資料
    mpc_data = pd.read_csv(file_path, sep='\s+', skiprows=43, header=None, usecols=[3, 5, 6],
                           names=["Epoch", "RA_deg", "Dec_deg"])
    # 篩選出合法的緯度範圍 (-90 到 90 度)
    mpc_data = mpc_data[(mpc_data["Dec_deg"] >= -90) & (mpc_data["Dec_deg"] <= 90)]
    return mpc_data

# 檔案中心的座標和時間
center_ra = 15 * 15  # 15h轉換為度數（15度/h）
center_dec = 0.5  # 直接使用度數
center_coord = SkyCoord(ra=center_ra * u.degree, dec=center_dec * u.degree, frame='icrs')

# 設定小行星資料的搜尋範圍
def find_asteroids_in_radius(mpc_data, center_coord, max_radius=5.0):
    for radius in range(1, int(max_radius) + 1):
        asteroid_coords = SkyCoord(ra=mpc_data["RA_deg"].values * u.degree, dec=mpc_data["Dec_deg"].values * u.degree)
        sep = center_coord.separation(asteroid_coords)
        nearby_asteroids = mpc_data[sep < radius * u.degree]
        
        if not nearby_asteroids.empty:
            print(f"在半徑 {radius} 度範圍內找到小行星")
            return nearby_asteroids, radius

    print("在給定的最大半徑內未找到小行星")
    return None, None

# 主程式
mpc_data = load_mpc_data('MPCORB.DAT')
asteroids_found, radius_used = find_asteroids_in_radius(mpc_data, center_coord)

if asteroids_found is not None:
    print("找到的小行星：")
    print(asteroids_found)
else:
    print("未找到小行星")
