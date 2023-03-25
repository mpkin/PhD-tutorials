# plotter.py: plots the interpolated solution, exact solution, and error

import numpy as np
import matplotlib.pyplot as plt

# Define list of data file names
data_files = ['interp.dat', 'exact.dat', 'error.dat']

# Loop over data files
for data_file in data_files:

    # Read data from file
    data = np.loadtxt(data_file)

    # Extract x, y, and z data
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    # Define grid for pcolormesh
    X, Y = np.meshgrid(np.unique(x), np.unique(y))
    Z = z.reshape(len(np.unique(y)), len(np.unique(x)))

    # Plot data
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.title(data_file)
    plt.pcolormesh(X, Y, Z, cmap='jet')
    plt.colorbar()
    plt.show()
