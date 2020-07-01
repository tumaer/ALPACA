#!/usr/bin/env python3

import os
import numpy as np
from scipy.signal import find_peaks
import h5py

def check_solution( resultsFolder, pythonScriptsFolder ) :

    resultsPath = os.path.join(resultsFolder, "OscillatingDrop/domain")

    files = []
    time = []
    for file in os.listdir(resultsPath):
        if file.endswith("h5"):
                files.append(file)
                time.append(float(file[5:13]))
    indices = np.argsort(time)
    time = np.array(time)[indices]
    files = np.array(files)[indices]

    total_energy = []
    for filename in files:

        with h5py.File(os.path.join(resultsPath, filename), "r") as data:
            velocityX = np.array(data["cell_data"]["velocity"][:,0,0])
            velocityY = np.array(data["cell_data"]["velocity"][:,1,0])
            levelset = np.array(data["cell_data"]["levelset"][:,0,0])

        energy = (velocityX**2 + velocityY**2)
        indices = np.where(levelset > 0)[0]
        energy = np.sum(energy[indices])
        total_energy.append(energy)

    # find the local maxima of kinetic energy, from which the oscillation period is computed
    total_energy = np.array(total_energy)

    indices, _ = find_peaks(total_energy, width=10)
    time_maxima = time[indices]

    T = np.diff(time_maxima)
    T = np.mean(T)

    R = 0.4
    rho_l = 100.0
    rho_g = 5.0
    sigma = 200.0
    oscillation_period = np.pi * np.sqrt((rho_l + rho_g) * np.power(R,3) / 6 / sigma)

    rel_error = (T - oscillation_period)/oscillation_period

    return abs(rel_error)
