# sez_to_ecef.py
#
# Usage: Converts SEZ vector components to ECEF
# Parameters:
#  o_lat_deg: Latitude in degrees
#  o_lon_deg: Longitude in degrees
#  o_hae_km:  Height Above Ellipsoid in km
#  s_km: South component in km
#  e_km: East component in km
#  z_km: Z component in km
# Output:
#  ECEF x, y & z-component in km
#
# Written by Tejas Vinod
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.
#
# Test Inputs:
# python3 sez_to_ecef.py 40.496 -80.246 0.37 0 1 0.3
# Expected Answer
# 822.957, -4787.242, 4120.457

# import Python modules
import math # math module
import sys # argv

# "constants"
R_E_KM = 6378.1363
e_E = 0.081819221456

# helper functions

## calc_denom
def calc_denom(ecc,lat_rad):
  return math.sqrt(1.0-ecc**2.0 * math.sin(lat_rad)**2.0)

def matrix_multiply(matrix1, matrix2):
    # Get the dimensions of the matrices
    rows_matrix1 = len(matrix1)
    cols_matrix1 = len(matrix1[0])
    rows_matrix2 = len(matrix2)
    cols_matrix2 = len(matrix2[0])
    
    # Check if matrices are able to be multiplied
    if cols_matrix1 != rows_matrix2:
        raise ValueError("Matrix multiplication is not possible. Columns of matrix1 must equal rows of matrix2.") 
    
    # Perform the multiplication
    result = [[0 for _ in range(cols_matrix2)] for _ in range(rows_matrix1)]
    for i in range(rows_matrix1):
        for j in range(cols_matrix2):
            for k in range(cols_matrix1):
                result[i][j] += matrix1[i][k] * matrix2[k][j]
    
    return result

def add_vectors(vector1, vector2):
    result = []
    for i in range(len(vector1)):
        # Check if elements are lists or just floats
        if isinstance(vector1[i], list):
            result.append(vector1[i][0] + vector2[i][0])
        else:
            result.append(vector1[i] + vector2[i][0])  # vector2 is likely still a list of lists
    
    return result

# initialize script arguments
o_lat_deg = float('nan') # Latitude in degrees
o_lon_deg = float('nan') # Longitude in degrees
o_hae_km  = float('nan') # Height Above Ellipsoid in km
s_km = float('nan') # South component in km
e_km = float('nan') # East component in km
z_km = float('nan') # Z component in km

# parse script arguments
if len(sys.argv)==7:
  o_lat_deg = float(sys.argv[1])
  o_lon_deg = float(sys.argv[2])
  o_hae_km  = float(sys.argv[3])
  s_km = float(sys.argv[4])
  e_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
   'Usage: '\
   'python3 o_lat_deg o_lon_deg o_hae_km s_km e_km z_km'\
  )
  exit()

# write script below this line
# Finding r_ecef of plane from origin
r_sez = [[s_km], [e_km], [z_km]]
o_lat_rad = o_lat_deg * math.pi/180
o_lon_rad = o_lon_deg * math.pi/180

Ry = [ 
    [ math.sin(o_lat_rad), 0, math.cos(o_lat_rad)],
    [                   0, 1,                   0],
    [-math.cos(o_lat_rad), 0, math.sin(o_lat_rad)] ] # ry(90-phi)

Rz = [ 
    [math.cos(o_lon_rad), -math.sin(o_lon_rad), 0],
    [math.sin(o_lon_rad),  math.cos(o_lon_rad), 0],
    [                  0,                    0, 1] ] # rz(theta)

r_ecef_1 = matrix_multiply(Rz,matrix_multiply(Ry,r_sez))

# Finding r_ecef of origin
r_ecef_o = [[0], [0], [0]]
denom = calc_denom(e_E,o_lat_rad)
C_E = R_E_KM/denom
S_E = R_E_KM*(1-e_E**2)/denom
r_ecef_o[0] = (C_E+o_hae_km)*math.cos(o_lat_rad)*math.cos(o_lon_rad)
r_ecef_o[1] = (C_E+o_hae_km)*math.cos(o_lat_rad)*math.sin(o_lon_rad)
r_ecef_o[2] = (S_E+o_hae_km)*math.sin(o_lat_rad)

r_ecef = add_vectors(r_ecef_o,r_ecef_1)

print('ecef_x_km: ' + str(r_ecef[0]))
print('ecef_y_km: ' + str(r_ecef[1]))
print('ecef_z_km: ' + str(r_ecef[2]))   