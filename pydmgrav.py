#/usr/bin/python
import time
from itertools import repeat, starmap
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import scipy.stats
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from tqdm import tqdm

# need to run ipcluster start
# before using this code!

# make blacklist
bl_path = Path("blacklist.txt")
bl = bl_path.read_text().split("\n")
BLACKLIST = [Path(i) for i in bl if (i != '') and (i[0] != "#")]


# make fig3 function
def make_fig3(path="deltag_result.txt"):
    fig3_t, fig3_y = np.loadtxt(path).T
    fig3_y = fig3_y * 1e9 * 1e5
   
    return interp1d(fig3_t,
                    fig3_y,
                    fill_value=0,
                    bounds_error=False,
                    kind="linear")


def load_file(path, minsize=100, dump_to_npy=False):
    """
    Loads an IGETS data product file at path
    """
    try:
        raw_data = np.genfromtxt(path,
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
    split_mask = np.isnan(raw_data[:, 0])\
                      + np.isnan(raw_data[:, 1])\
                      + np.isnan(raw_data[:, 2])\
                      + np.isclose(raw_data[:, 0], 77777777)\
                      + np.isclose(raw_data[:, 0], 88888888)\
                      + np.isclose(raw_data[:, 0], 99999999)
    splits = np.where(split_mask)[0]

    datalist = []
    if splits.size != 0:
        print("found a split!")
        print(path)
        print(splits)
        for count, index in enumerate(splits):
            if count == 0:
                if index > minsize:
                    datalist.append(raw_data[:index])
            else:
                # split_mask[count - 1] + 1 is the index of the next data point
                # after the last bad one
                candidate = raw_data[splits[count - 1] + 1:index]
                if len(candidate) > minsize:
                    datalist.append(candidate)
    else:
        datalist = [raw_data]

    # pull out just the year/month/day column and the hour/minute/second column
    # and cast them to strings (int as an intermediate to remove spurious
    # decimals
    resultlist = []
    for data in datalist:
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
        resultlist.append(result)
    if dump_to_npy:
        # strip the extension
        new_path = path[:-4]
        for i, r in enumerate(resultlist):
            np.save(new_path + "_" + str(i) + ".npy", r)
        return

    fft_list = np.array(list(map(do_fft_on_data, resultlist)))

    mean_fft = fft_list.mean(axis=0)
    if type(mean_fft) is not np.ndarray:
        return ("issues with", path)
    else:
        return mean_fft


def gaussian(f, a, sigma, f0):
    return a * np.exp(-(((f - f0) / sigma)**2))


def lorentzian(f, A, gamma, f0):
    return (A / np.pi) * (gamma / 2) / ((f - f0)**2 + (gamma / 2)**2)


def split_lor(f, A, gamma, f0, m, b):
    split = 1 / (60 * 60 * 24)
    return lorentzian(f, A / 2, gamma, f0 - split) + lorentzian(
        f, A / 2, gamma, f0 + split) + m * f + b


# the figure 3 function
fig3_f = make_fig3()


def one_on_f(f, A, B, C, D, y0):
    return A / f + B / f**2 + C / f**3 + D / f**4 + y0


def do_fft_on_data(data,
                   max_freq=1 / 180,
                   min_freq=1 / 172800,
                   bl_tidal_cutoff=3.655e-5,
                   interp_freq_step=1.87e-07,
                   inject_amplitude=None):
    if type(data) is not np.ndarray:
        data = np.load(data)

    if inject_amplitude is not None:
        data[:, 1] = data[:, 1] + inject_amplitude * fig3_f(data[:, 0] - data[0, 0]) / 2

    timestep = data[1, 0] - data[0, 0]
    freqs = np.fft.rfftfreq(len(data), timestep)
    fft = np.fft.rfft(data[:, 1])
    psd = (timestep**2) * (np.abs(fft)**2) / (data[-1, 0] - data[0, 0])
    # interpolate the spectra so we can average them properly later
    interp_spectra = interp1d(freqs, psd)

    # the nyquist frequency
    # max_freq = 1 / (3 * 60)
    # min freq corresponds to a 200 min period
    # min_freq = 1 / (200 * 60)
    # frequency step, use 2x the highest frequency step ive seen
    # to avoid aliasing
    # interp_freq_step=1.87e-07
    new_freqs = np.arange(min_freq, max_freq, interp_freq_step)

    new_psd = interp_spectra(new_freqs)
    return new_psd


def convert_files(level=3, basedir="./"):
    """
    level can be 2, 3, or "both"
    """
    basedir = Path(basedir)
    files = []
    files3 = list(basedir.rglob("Level3/**/*RESMIN*.ggp"))
    files2 = list(basedir.rglob("Level2/**/*CORMIN*.ggp"))

    if level == 3:
        files += files3
    if level == 2:
        files += files2
    if level == "both":
        files = files3 + files2

    for f in tqdm(files):
        x = load_file(f, dump_to_npy=True)
    return


def main_npy(level=3,
             basedir="./",
             subtract_global_baseline=True,
             max_freq=1 / 180,
             min_freq=1 / 172800,
             bl_tidal_cutoff=3.655e-5,
             interp_freq_step=1.87e-07,
             parallel=True,
             inject_amplitude=None):

    basedir = Path(basedir)
    tstart = time.time()
    if level == 3:
        files = basedir.rglob("Level3/**/*RESMIN*.npy")
    elif level == 2:
        files = basedir.rglob("Level2/**/*CORMIN*.npy")
    files = [i for i in files if i not in BLACKLIST]

    args = zip(files, repeat(max_freq), repeat(min_freq),
               repeat(bl_tidal_cutoff), repeat(interp_freq_step),
               repeat(inject_amplitude))
    # load the files in parallel
    print("start loading files")
    if parallel:
        with Pool(processes=8) as p:
            psds = list(
                tqdm(p.starmap(do_fft_on_data, args, chunksize=50),
                     total=len(files)))
    else:
        psds = list(starmap(do_fft_on_data, tqdm(files)))

    psds = np.array(psds)
    # single threaded version of above for debugging:
    # ffts = list(map(load_file, files))
    print("done loading files", time.time() - tstart)

    # generate the frequency array
    freqs = np.arange(min_freq, max_freq, interp_freq_step)
    mean_psd = np.mean(psds, axis=0)
    mean_psd_sigma = scipy.stats.sem(psds, axis=0)
    mean_asd = np.sqrt(mean_psd)
    # error propigation for square root
    mean_asd_sigma = (1 / 2) * (mean_psd_sigma / mean_psd) * mean_asd

    if subtract_global_baseline:
        # remove the baseline
        # only fit data that's higher in freqeuncy than the tides
        tide_mask = freqs > bl_tidal_cutoff
        fit_results = curve_fit(one_on_f,
                                freqs[tide_mask],
                                mean_asd[tide_mask],
                                sigma=mean_psd_sigma[tide_mask],
                                absolute_sigma=True)
        # subtract the baseline
        mean_asd = mean_asd - one_on_f(freqs, *fit_results[0])

    print("done averaging", time.time() - tstart)
    return np.array((freqs, mean_asd, mean_asd_sigma))


def main_raw(level=3,
             basedir="./",
             subtract_global_baseline=True,
             max_freq=1 / 180,
             min_freq=1 / 172800,
             bl_tidal_cutoff=3.655e-5,
             interp_freq_step=1.87e-07):
    basedir = Path(basedir)
    tstart = time.time()
    if level == 3:
        files = basedir.rglob("Level3/**/*RESMIN*.ggp")
    elif level == 2:
        files = basedir.rglob("Level2/**/*CORMIN*.ggp")
    files = [i for i in files if i not in BLACKLIST]
    # load the files in parallel
    print("start loading files")
    with Pool(processes=8) as p:
        psds = p.map(load_file, files, chunksize=5)
    # single threaded version of above for debugging:
    # ffts = list(map(load_file, files))
    print("done loading files", time.time() - tstart)

    # Print and then remove any failed files

    print("Failed Files:")
    failed_files = [i for i in psds if type(i[0]) is str]
    for i in failed_files:
        print(i)
    # remove the failed files
    psds[:] = [i for i in psds if type(i[0]) is not str]
    [print(i) for i in psds if type(i[0]) is str]
    psds = np.array([i for i in psds if i is not None])
    # generate the frequency array
    freqs = np.arange(min_freq, max_freq, interp_freq_step)

    mean_psd = np.mean(psds, axis=0)
    mean_psd_sigma = scipy.stats.sem(psds, axis=0)
    mean_asd = np.sqrt(mean_psd)
    # error propigation for square root
    mean_asd_sigma = (1 / 2) * (mean_psd_sigma / mean_psd) * mean_asd

    if subtract_global_baseline:
        # remove the baseline
        # only fit data that's higher in freqeuncy than the tides
        tide_mask = freqs > bl_tidal_cutoff
        fit_results = curve_fit(one_on_f,
                                freqs[tide_mask],
                                mean_asd[tide_mask],
                                sigma=mean_psd_sigma[tide_mask],
                                absolute_sigma=True)
        # subtract the baseline
        mean_asd = mean_asd - one_on_f(freqs, *fit_results[0])

    print("done averaging", time.time() - tstart)
    return np.array((freqs, mean_asd, mean_asd_sigma))


def do_all_sites():
    """
    This function does the data reduction on each site, one at a time,
    then plots them and returns the results in an array.
    """
    p = Path(".")
    sites = [x for x in p.iterdir() if x.is_dir() and str(x)[0].isupper()]

    # filter out sites without level3 data
    sites = [i for i in sites if list(i.rglob("/Level3"))]

    results = []
    for i, site in enumerate(sites):
        print(site)
        result = main_npy(basedir=site)
        plt.plot(result[0], result[1], label=site)
        plt.legend()
        plt.pause(.1)
        results.append(result)
    return np.array(results)


def dig_single_site(basedir):
    """
    Crappy little function to dig into specific sites and look at the
    individual 1-month spectra. Mostly a scratchpad function, as what needs
    investigating varies.
    """
    basedir = Path(basedir)
    files = basedir.rglob("Level3/**/*RESMIN*.npy")

    # the nyquist frequency
    max_freq = 1 / (3 * 60)
    # min freq corresponds to a 200 min period
    min_freq = 1 / (48 * 60 * 60)
    # frequency step, use 2x the highest frequency step ive seen
    # to avoid aliasing
    interp_freq_step = 1.87e-07
    freqs = np.arange(min_freq, max_freq, interp_freq_step)
    max_val = 0

    mask = (freqs > 0.003) * (freqs < 0.004)

    for f in files:
        data = np.load(f)
        psd = do_fft_on_data(data)
        plt.plot(freqs, psd)
        if np.std(psd[mask]) > max_val:
            max_val = np.std(psd[mask])
        print(f)
        print(np.std(psd[mask]))


def rchisqr(r, p, mask):
    residuals = split_lor(r[0][mask], *p) - r[1][mask]
    chisqr = np.sum((residuals / r[2][mask]) ** 2)
    return chisqr / (len(r[0]) - len(p))


def do_injection():
    inject_amplitudes = np.logspace(np.log(.001), np.log(10), 20, base=np.e)
    rs = [main_npy(inject_amplitude=i) for i in inject_amplitudes]
    p0 = [7.31462455e-05,  5.09127562e-07,  3.03933293e-04, -8.39103118e+04, 3.92960778e+01]
    ps = []
    rcs = []
    for r in rs:
        mask = (r[0] > .000275) * (r[0] < .000325)
        p, cov = curve_fit(split_lor, r[0][mask], r[1][mask], p0=p0, sigma=r[2][mask])
        ps.append(p)
        rcs.append(rchisqr(r, p, mask))
    import ipdb; ipdb.set_trace()
    return (ps, rcs)


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
