#!/usr/bin/env python
"""
This module provides functionality to read data from 'data.temp' and 'data.spec' files
in a specified directory and plot the data using the Lomb-Scargle periodogram.

Functions:
- read_file(cwd: str, filename: str) -> tuple[np.ndarray, np.ndarray]: Reads data from a
    file in a specified directory and returns two arrays.
- plot(time, amplitude, freq, power): Plots the data using the Lomb-Scargle periodogram.

Usage:
The main() function of this module loads data from 'data.temp' and 'data.spec' files in
a specified directory and plots the data using the Lomb-Scargle periodogram.

Example:
    To run the main() function, import the module and call it:
        import lomb_scargle_periodogram
        lomb_scargle_periodogram.main()
"""

import subprocess
import os

import matplotlib.pyplot as plt
import numpy as np
import colorama

def read_file(cwd: str, filename: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Reads a CSV file with a header row and two columns of numerical data separated by semicolons,
    located in the current working directory specified by `cwd` and `filename`, and returns
    the data as two NumPy arrays.

    Parameters:
    cwd : str
        The current working directory where the CSV file is located.
    filename : str
        The name of the CSV file to be read.

    Returns:
    Tuple[np.ndarray, np.ndarray]
        Two 1D arrays containing the data from the columns of the file.

    Raises:
    FileNotFoundError
        If the specified directory or file does not exist.
    ValueError
        If the number of entries in the columns of the file do not match.
    """
    if not os.path.isdir(cwd):
        raise FileNotFoundError(f"The specified directory '{cwd}' does not exist.")

    file_path = os.path.join(cwd, filename)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The specified file '{filename}' does not exist in '{cwd}' directory.")

    col_1, col_2 = np.loadtxt(file_path, dtype=float, comments='#', delimiter=';', skiprows=1, unpack=True)

    # sanity check
    if len(col_1) != len(col_2):
        raise ValueError(f"The number of entries in columns of file '{filename}' do not match.")

    return col_1, col_2

def plot(time: np.ndarray, amplitude:np.ndarray, freq: np.ndarray, power: np.ndarray) -> None:
    """
    Plots a Lomb-Scargle periodogram, given the time and amplitude arrays of a signal and the corresponding
    frequency and power arrays obtained from the Lomb-Scargle periodogram analysis.

    Parameters:
    time : numpy.ndarray
        The time values of the signal, in seconds.
    amplitude : numpy.ndarray
        The amplitude values of the signal.
    freq : numpy.ndarray
        The frequency values obtained from the Lomb-Scargle periodogram analysis, in Hz.
    power : numpy.ndarray
        The power values obtained from the Lomb-Scargle periodogram analysis.

    Returns:
    None
    """
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    fig.canvas.manager.set_window_title('Lomb-Scargle periodogram')

    ax1.set_xlabel(r'time [$s$]')
    ax1.set_ylabel(r'amplitude [$arb.u.$]')
    ax2.set_xlabel(r'frequency [$1/s$]')
    ax2.set_ylabel(r'power [$arb.u.$]')

    ax1.plot(time, amplitude, '.', color='Cornflowerblue')
    ax2.plot(freq, power, color='Cornflowerblue')

    ax1.set_xlim([time[0]-1, time[-1]+2])
    ax2.set_xlim([freq[0]-.003, freq[-1]])

    # significance levels
    plt.axhline(14.51, 0, .95, color='Maroon', linewidth=.3)  # p = 0.0001
    plt.axhline(12.21, 0, .95, color='Maroon', linewidth=.3)  # p = 0.001
    plt.axhline(09.90, 0, .95, color='Maroon', linewidth=.3)  # p = 0.01
    plt.axhline(07.55, 0, .95, color='Maroon', linewidth=.3)  # p = 0.1

    # significance labels
    ax2.text(freq[-1]*.96, 14.51, '.0001', size=5)
    ax2.text(freq[-1]*.96, 12.21, '.001', size=5)
    ax2.text(freq[-1]*.96,  9.90, '.01', size=5)
    ax2.text(freq[-1]*.96,  7.55, '.1', size=5)

    plt.tight_layout()
    plt.show()

def main() -> None:
    """
    Load data from 'data.temp' and 'data.spec' files in the specified directory.
    Plot the data using the Lomb-Scargle periodogram.
    """
    colorama.init(autoreset=True)
    red = colorama.Fore.RED + colorama.Style.BRIGHT

    cmd = 'pwd'
    cwd = subprocess.run(cmd, shell=True, encoding='utf-8', stdout=subprocess.PIPE, check=True).stdout.strip()

    # Load data from 'data.temp' and 'data.spec' files
    try:
        time, amplitude = read_file(cwd, 'data.temp')
        freq, power = read_file(cwd, 'data.spec')
    except FileNotFoundError as e:
        print(red + str(e))
        return
    except ValueError as e:
        print(red + str(e))
        return

    # Plot the data using the Lomb-Scargle periodogram
    plot(time, amplitude, freq, power)

if __name__ == '__main__':
    main()
