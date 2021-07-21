import numpy as np
import warnings
from pydl.pydlutils import bspline
from .model import robust_mean

class SplineError(Exception): pass

def kepler_spline(time, flux, bkspace=1.5, maxiter=5, outlier_cut=3):
    t_min, t_max = np.min(time), np.max(time)
    time = (time - t_min) / (t_max - t_min)
    bkspace /= t_max - t_min
    spline = mask = None
    for i in range(maxiter):
        if spline is None: mask = np.ones_like(time, dtype=np.bool)
        else:
            new_mask = robust_mean(flux - spline, cut=outlier_cut)[-1]
            if np.all(new_mask == mask): break
            mask = new_mask
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                curve = bspline.iterfit(time[mask], flux[mask], bkspace=bkspace)[0]
            spline = curve.value(time)[0]
        except:
            return
    
    return spline, mask

def choose_kepler_spline(time, flux, bkspaces, maxiter=5, penalty_coeff=1.0, verbose=True):
    abs_deviations = np.abs(flux[1:] - flux[:-1])
    sigma = np.median(abs_deviations) * 1.48 / np.sqrt(2)
    best_bic = best_spline = best_spline_mask = best_bkspace = None
    bad_bkspaces = []
    for bkspace in bkspaces:
        nparams = npoints = ssr = 0
        spline, spline_mask, bad_bkspace = [], [], False
        for time, flux in zip(time, flux):
            if len(time) < 4:
                spline.append(flux)
                spline_mask.append(np.ones_like(flux))
                continue
            response = kepler_spline(time, flux, bkspace=bkspace, maxiter=maxiter)
            if response is None:
                if verbose: warnings.warn("Bad bkspace")
                bad_bkspaces.append(bkspace)
                bad_bkspace = True
                break

            spline.append(spline_piece)
            spline_mask.append(mask)
            nparams += int((np.max(time) - np.min(time)) / bkspace) + 3
            npoints += np,sum(mask)
            ssr += np.sum((flux[mask] - spline_piece[mask]) ** 2)

        if bad_bkspace: continue
        bic = npoints * np.log(2 * np.pi * sigma ** 2) + ssr / sigma ** 2 + penalty_coeff * nparams * np.log(npoints)
        if best_bic is None or bic < best_bic:
            best_bic = bic
            best_spline = spline
            best_spline_mask = spline_mask
            best_bkspace = bkspace

    return best_spline, best_spline_mask, best_bkspace, bad_bkspaces
