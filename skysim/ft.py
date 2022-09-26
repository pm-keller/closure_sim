"""

ft.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 14/02/2022

Description:
    Fourier Transform Tools

"""

import numpy as np
from scipy.signal import windows, csd


def fft(v, W=None, norm=None, inverse=False, shift=False, axis=-1):
    """Compute FFT

    Parameters
    ----------
    v : ndarray
        input data
    W : ndarray, optional
        window function, by default None
    norm : str, optional
        normalisation (None or "ortho"), by default None
    inverse : bool, optional
        If True, perform inverse FFT, by default False
    shift : bool, optional
        If True, shift FFT to centre, by default False
    axis : int, optional
        axis along which to compute fourier transform, by default -1

    Returns
    -------
    ndarray
        Fourier transformed data
    """

    v = np.moveaxis(v, axis, -1)

    # generate boxcar window
    if W is None:
        W = np.ones(v.shape[-1])

    # compute inverse FFT
    if inverse:
        ft = np.fft.ifft(W * v, axis=-1, norm=norm)

        # shift FFT to centre
        if shift:
            ft = np.fft.ifftshift(ft, axes=-1)

    # compute FFT
    else:
        ft = np.fft.fft(W * v, axis=-1, norm=norm)

        # shift FFT to centre
        if shift:
            ft = np.fft.fftshift(ft, axes=-1)

    return np.moveaxis(ft, -1, axis)


def xps(v1, v2, fs=1.0, W=None, inverse=False, shift=False, axis=-1, **kwargs):
    """Compute FFT Cross-Power Spectrum

    Parameters
    ----------
    v1 : ndarray
        input data 1
    v2 : ndarray
        input data 2
    fs : float, optional
        sampling frequency, by default 1.0
    W : ndarray, optional
        window function, by default None
    inverse : bool, optional
        If True, perform inverse FFT, by default False
    shift : bool, optional
        If True, shift FFT to centre, by default False
    axis : int, optional
        axis along which to compute fourier transform, by default -1

    Returns
    -------
    ndarray
        Cross power spectrum between input data 1 and 2
    """

    # compute FFT
    ft1 = fft(v1, W, inverse=inverse, shift=shift, axis=axis, **kwargs)
    ft2 = fft(v2, W, inverse=inverse, shift=shift, axis=axis, **kwargs)

    # inverse normalisation
    if inverse:
        ft1 = ft1 * ft1.shape[-1]
        ft2 = ft2 * ft2.shape[-1]

    return ft1 * ft2.conjugate() / fs ** 2
