#/usr/bin/python
import glob
import time

import ipyparallel as ipp
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# need to run ipcluster start
# before using this code!


def load_file(path):
    """
    Loads an IGETS data product file at path
    """
    try:
        data = np.genfromtxt(path,
                             skip_header=30,
                             skip_footer=2,
                             invalid_raise=False,
                             usecols=(0, 1, 2),
                             encoding='latin1',
                             delimiter=(9, 6, 10))
    except UnicodeDecodeError:
        return ("Unicode Issue", path)
        # raise ValueError("path is  " + path)

    # pathalogical lines exist two ways, either they start with a nan
    # or with a number like 77777777 or 99999999. filter both
    #
    data = data[~np.isnan(data[:, 0])]
    data = data[~np.isclose(data[:, 0], 77777777)]
    data = data[~np.isclose(data[:, 0], 88888888)]
    data = data[~np.isclose(data[:, 0], 99999999)]
    # pull out just the year/month/day column and the hour/minute/second column
    # and cast them to strings (int as an intermediate to remove spurious
    # decimals
    yyyymmdd = data[:, 0].astype(int).astype(str)
    hhmmss = data[:, 1].astype(int).astype(str)
    hhmmss = np.core.defchararray.zfill(hhmmss, 6)

    years = [i[0:4] for i in yyyymmdd]
    months = [i[4:6] for i in yyyymmdd]
    days = [i[6:8] for i in yyyymmdd]

    hours = [i[0:2] for i in hhmmss]
    minutes = [i[2:4] for i in hhmmss]
    seconds = [i[4:6] for i in hhmmss]

    ymd = ["-".join(i) for i in zip(years, months, days)]
    hms = [":".join(i) for i in zip(hours, minutes, seconds)]
    try:
        timestamps = np.array(["T".join(i) for i in zip(ymd, hms)],
                              dtype=np.datetime64)
    except ValueError:
        print("something wrong with date.")
        print(path)
        # raise ValueError("path is: " + path)
        print("Skipping this file...")
        return ("date issue with", path)

    result = np.column_stack((timestamps.astype(np.float), data[:, 2]))
    return do_fft_on_data(result)


def do_fft_on_data(data, max_freq=1 / 180, min_freq=1 / 172800):
    timestep = data[1, 0] - data[0, 0]
    freqs = np.fft.rfftfreq(len(data), timestep)
    fft = np.abs(np.fft.rfft(data[:, 1]))

    # remove the baseline
    # below this frequency we run into the tides
    bl_tidal_cutoff = 3.655e-5

    one_on_f = lambda f, A, B, C, D, y0: A / f + B / f**2 + C / f**3 + D / f ** 4 +  y0
    tide_mask = freqs > bl_tidal_cutoff
    # plt.plot(freqs[tide_mask], fft[tide_mask])
    # fit only the non-tidal region
    fit_results = curve_fit(one_on_f,
                            freqs[tide_mask],
                            fft[tide_mask])
    print(fit_results[0])
    # subtract the baseline
    fft = fft - one_on_f(freqs, *fit_results[0])
    # plt.plot(freqs[tide_mask], fft[tide_mask])
    # plt.pause(.5)
    # plt.clf()

    interp_spectra = interp1d(freqs, np.abs(fft))

    # the nyquist frequency
    #max_freq = 1 / (3 * 60)
    # min freq corresponds to a 200 min period
    #min_freq = 1 / (200 * 60)
    # frequency step, use 2x the highest frequency step ive seen
    # to avoid aliasing
    interp_freq_step = 1.87e-07
    new_freqs = np.arange(min_freq, max_freq, interp_freq_step)

    new_fft = interp_spectra(new_freqs)
    return new_fft


def main(level=3, basedir="./"):
    tstart = time.time()
    if level == 3:
        files = glob.glob(basedir + "**/Level3/**/*RESMIN*.ggp",
                          recursive=True)
    elif level == 2:
        files = glob.glob(basedir + "**/Level2/**/*CORMIN*.ggp",
                          recursive=True)
    # interpolate results to add them
    # load the files in parallel
    print("start loading files")
    rc = ipp.Client()
    dview = rc[:]

    ffts = dview.map_sync(load_file, files)
    # ffts = list(map(load_file, files))
    print("done loading files", time.time() - tstart)

    # remove any failed files

    print("Failed Files:")
    failed_files = [i for i in ffts if type(i[0]) is str]
    for i in failed_files:
        print(i)
    # remove the failed files
    ffts[:] = [i for i in ffts if type(i[0]) is not str]

    [print(i) for i in ffts if type(i[0]) is str]
    ffts = np.array([i for i in ffts if i is not None])
    #import ipdb; ipdb.set_trace()
    # the nyquist frequency
    max_freq = 1 / (3 * 60)
    # min freq corresponds to a 200 min period
    min_freq = 1 / (172800)
    # frequency step, use 2x the highest frequency step ive seen
    # to avoid aliasing
    interp_freq_step = 1.87e-07
    freqs = np.arange(min_freq, max_freq, interp_freq_step)

    mean_fft = np.mean(ffts, axis=0)
    print("done averaging", time.time() - tstart)
    return freqs, mean_fft


def plot_results(freqs, fft):
    periods = 1 / (60 * freqs)

    plt.subplot(211)
    plt.plot(freqs, fft)
    plt.xlabel("F (Hz)")
    plt.ylabel("Acceleration / Hz (in nm/s)")
    plt.subplot(212)
    plt.plot(periods, fft)
    plt.ylabel("Acceleration / Hz (in nm/s)")
    plt.xlabel("Period (in Minutes)")
