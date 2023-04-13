import os, sys
import datetime
from datetime import datetime
from datetime import timedelta
import numpy as np

from mpi4py import MPI
import shutil
import scipy.stats
import matplotlib.pyplot as plt

var_names = ["precipitation_bilinear", "precipitation_conserve", "temperature",
             "solar_radiation", "longwave_radiation",
             "specific_humidity", "wind_speed", "surface_pressure"]
var_units = ["mm/s", "mm/s", "K", "W/m2", "W/m2", "kg/kg", "m/s", "Pa"]
# only for precip, temp, SWRad, LWRad
std_dev = [0.5, 0.3, 50.0]     # 2.0,[0.2, 0.2, 2, 50, 80, 0.1, 0.1, 100]
mean_var= [1.0, 1.0, 0.0]      # 0.0,
multi_bool = [True, True, False]   # False, [True, True, False, False, False, True, True, False]
# only for          precip, SWRad, LWRad, temp,
pert_corr_matrix = [[1.0, -0.8, 0.5],
                    [-0.8, 1.0, -0.5],
                    [ 0.5, -0.5, 1.0]]
temp_corr_hr = 24.0
spat_corr_km = 50.0
ens_size = 20

noise = scipy.stats.multivariate_normal.rvs(mean=np.zeros(3), cov=pert_corr_matrix, size=ens_size)
print(noise)
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# ax2.contourf(noise.reshape((4,4)))
# plt.show()


#### Use convolution of random noise with gausian kernel--faster than correlated noise generation
# https://stackoverflow.com/questions/63816481/faster-method-for-creating-spatially-correlated-noise
# import numpy as np
# import scipy.signal
# import matplotlib.pyplot as plt
#
# # Compute filter kernel with radius correlation_scale (can probably be a bit smaller)
# correlation_scale = 50
# x = np.arange(-correlation_scale, correlation_scale)
# y = np.arange(-correlation_scale, correlation_scale)
# X, Y = np.meshgrid(x, y)
# dist = np.sqrt(X*X + Y*Y)
# filter_kernel = np.exp(-dist**2/(2*correlation_scale))
#
# # Generate n-by-n grid of spatially correlated noise
# n = 50
# noise = np.random.randn(n, n)
# noise = scipy.signal.fftconvolve(noise, filter_kernel, mode='same')
# plt.contourf(np.arange(n), np.arange(n), noise)
# plt.savefig("fast.png")
#
# ### slow one
#
# import scipy.stats
# import numpy as np
# import scipy.spatial.distance
# import matplotlib.pyplot as plt
#
# # Create a 50-by-50 grid; My actual grid will be a LOT larger
# X,Y = np.meshgrid(np.arange(50),np.arange(50))
#
# # Create a vector of cells
# XY = np.column_stack((np.ndarray.flatten(X),np.ndarray.flatten(Y)))
#
# # Calculate a matrix of distances between the cells
# dist = scipy.spatial.distance.pdist(XY)
# dist = scipy.spatial.distance.squareform(dist)
#
# # Convert the distance matrix into a covariance matrix
# correlation_scale = 50
# cov = np.exp(-dist**2/(2*correlation_scale)) # This will do as a covariance matrix
#
# # Sample some noise !slow!
# noise = scipy.stats.multivariate_normal.rvs(
#         mean = np.zeros(50**2),
#         cov = cov)
#
# # Plot the result
# plt.contourf(X,Y,noise.reshape((50,50)))