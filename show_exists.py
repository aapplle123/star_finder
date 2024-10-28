import pandas as pd
from skyfield.api import load
from skyfield.data.mpc import load_mpcorb_dataframe

# List of asteroid data (with approximate RA and Dec)
asteroid_data = [
    {"name": "00001", "ra": 80.25414, "dec": 10.5879, "obs_date": "20240806"},
    {"name": "00002", "ra": 172.90614, "dec": 34.92186, "obs_date": "20240817"},
    {"name": "00003", "ra": 169.83829, "dec": 12.98815, "obs_date": "20240607"},
    {"name": "00004", "ra": 103.70471, "dec": 7.14399, "obs_date": "20240306"},
    {"name": "00005", "ra": 141.45931, "dec": 5.35896, "obs_date": "20240603"},
    {"name": "00006", "ra": 138.61983, "dec": 14.73436, "obs_date": "20240817"},
    {"name": "00007", "ra": 259.49897, "dec": 5.51948, "obs_date": "20240901"},
    {"name": "00008", "ra": 110.84574, "dec": 5.89031, "obs_date": "20240228"},
    {"name": "00009", "ra": 68.86867, "dec": 5.5782, "obs_date": "20240603"},
    {"name": "00010", "ra": 283.14844, "dec": 3.83162, "obs_date": "20240808"}
]

# Load timescale
ts = load.timescale()

def parse_date(date_str):
    """Convert the date string (YYYYMMDD) to a Skyfield time object."""
    year = int(date_str[:4])
    month = int(date_str[4:6])
    day = int(date_str[6:])
    return ts.utc(year, month, day)

def count_matching_asteroids(asteroid_data, mpcorb_data):
    matching_count = 0

    for asteroid in asteroid_data:
        asteroid_name = asteroid['name']
        asteroid_ra = asteroid['ra']
        asteroid_dec = asteroid['dec']
        asteroid_obs_date = asteroid['obs_date']

        print(f"Checking asteroid {asteroid_name} at RA: {asteroid_ra}, Dec: {asteroid_dec}, Date: {asteroid_obs_date}")

        # Filter the MPCORB data to check if any objects match the RA and Dec
        for index, row in mpcorb_data.iterrows():
            try:
                # Extract relevant orbital parameters
                inclination_mpc = float(row['Inclination'])  # Example: using Inclination as a substitute for Dec
                semi_major_axis_mpc = float(row['a'])  # Semi-major axis as an additional parameter
                obs_date_mpc = row['Reference'][-8:]  # Extract the date from 'Reference'

                # Check if the inclination and semi-major axis match approximately, and the observation date matches
                if abs(inclination_mpc - asteroid_dec) < 1.0 and obs_date_mpc == asteroid_obs_date:
                    print(f"Match found for asteroid {asteroid_name} in MPCORB data")
                    matching_count += 1
                    break  # Stop checking once a match is found
            except Exception as e:
                print(f"Error processing data: {e}")

    return matching_count

def main(mpcorb_file):
    # Load the minor planet data
    try:
        with open(mpcorb_file, 'rb') as f:
            mpcorb_data = load_mpcorb_dataframe(f)
            print("MPCORB data loaded successfully.")

            # Compare the asteroid data with the MPCORB data
            match_count = count_matching_asteroids(asteroid_data, mpcorb_data)
            print(f"\nTotal matching asteroids: {match_count}")

    except Exception as e:
        print(f"Error loading MPCORB data: {e}")

if __name__ == "__main__":
    mpcorb_file = "MPCORB.DAT"  # Update with the actual path to the MPCORB file
    main(mpcorb_file)
