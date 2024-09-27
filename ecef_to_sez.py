# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py ecef_x_km ecef_y_km ecef_z_km o_lat_deg o_lon_deg o_hae_km
# Converts ECEF vector components to SEZ
#
# Parameters:
# ecef_x_km ECEF x coordinate in km
# ecef_y_km ECEF y coordinate in km
# ecef_z_km ECEF z coordinate in km
# o_lat_deg Observer's latitude in deg
# o_lon_deg Observer's longitude in deg
# o_hae_km Observer's height above ellipsoid in km
# Output:
# Prints the SEZ coordinates (s_km, e_km, z_km)

# Written By: Samuel Jacobson

import math
import sys

# Constants
R_E_KM = 6378.137  # Earth's equatorial radius in km (consistent with your value)
E_E = 0.081819221456  # Earth's eccentricity

# Parse script arguments
if len(sys.argv) == 7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km = float(sys.argv[3])
    x_km = float(sys.argv[4])
    y_km = float(sys.argv[5])
    z_km = float(sys.argv[6])
else:
    print('Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km')
    exit()

# Compute delta_x, delta_y, delta_z
delta_x = x_km - o_x_km
delta_y = y_km - o_y_km
delta_z = z_km - o_z_km


# Function to calculate denominator
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0 - (ecc ** 2) * (math.sin(lat_rad) ** 2))


# Function to compute geodetic latitude, longitude, and height from ECEF coordinates
def ECEF_to_LLH(r_x_km, r_y_km, r_z_km):
    # Calculate longitude
    lon_rad = math.atan2(r_y_km, r_x_km)

    # Initialize latitude
    lat_rad = math.asin(r_z_km / math.sqrt(r_x_km ** 2 + r_y_km ** 2 + r_z_km ** 2))
    r_lon_km = math.sqrt(r_x_km ** 2 + r_y_km ** 2)
    prev_lat_rad = float('nan')

    # Iteratively find latitude
    c_E = float('nan')
    count = 0
    while (math.isnan(prev_lat_rad) or abs(lat_rad - prev_lat_rad) > 1e-7) and count < 5:
        denom = calc_denom(E_E, lat_rad)
        c_E = R_E_KM / denom
        prev_lat_rad = lat_rad
        lat_rad = math.atan((r_z_km + c_E * (E_E ** 2) * math.sin(lat_rad)) / r_lon_km)
        count += 1

    # Calculate height above ellipsoid (hae)
    hae_km = r_lon_km / math.cos(lat_rad) - c_E
    return lat_rad, lon_rad, hae_km


# Compute observer's geodetic latitude and longitude
lat_rad, lon_rad, o_hae_km = ECEF_to_LLH(o_x_km, o_y_km, o_z_km)

# Compute sin and cos of latitude and longitude
sin_lat = math.sin(lat_rad)
cos_lat = math.cos(lat_rad)
sin_lon = math.sin(lon_rad)
cos_lon = math.cos(lon_rad)

# Compute s, e, z components using the rotation matrix
s_km = (sin_lat * cos_lon) * delta_x + (sin_lat * sin_lon) * delta_y - (cos_lat) * delta_z
e_km = -sin_lon * delta_x + cos_lon * delta_y
z_km = (cos_lat * cos_lon) * delta_x + (cos_lat * sin_lon) * delta_y + (sin_lat) * delta_z

# Print results
print(s_km)
print(e_km)
print(z_km)