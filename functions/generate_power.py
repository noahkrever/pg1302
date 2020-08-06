def generate_power_DRW(N, dt, sigma, tau, generate_complex=False, random_state=None):
    import matplotlib.pyplot as plt
    from astroML.utils import check_random_state
    import numpy as np
    """Generate a power-law light curve
        This uses the method from Timmer & Koenig [1]_
           
        Parameters
           
        ----------
           
        N : integer
           
        Number of equal-spaced time steps to generate
           
        dt : float
           
        Spacing between time-steps
           
        beta : float
           
        Power-law index.  The spectrum will be (1 / f)^beta
           
        generate_complex : boolean (optional)
           
        if True, generate a complex time series rather than a real time series
           
        random_state : None, int, or np.random.RandomState instance (optional)
           
        random seed or random number generator
           
           
           
        Returns
           
        -------
           
        x : ndarray
           
        the length-N
           
           
           
        References
           
        ----------
           
        .. [1] Timmer, J. & Koenig, M. On Generating Power Law Noise. A&A 300:707
           
    """


    random_state = check_random_state(random_state)
    
    dt = float(dt)
    
    N = int(N)
    
    
    
    Npos = int(N / 2)
    
    Nneg = int((N - 1) / 2)
    
    domega = (2 * np.pi / dt / N)
    
    df = (1 / dt / N)
    
    
    
    if generate_complex:
        
        omega = domega * np.fft.ifftshift(np.arange(N) - int(N / 2))

    else:
        
        omega = domega * np.arange(Npos + 1)

    if generate_complex:
    
        freq = df * np.fft.ifftshift(np.arange(N) - int(N / 2))
    
    else:
        
        freq = df * np.arange(Npos + 1)


    x_fft = np.zeros(len(omega), dtype=complex)

    x_fft.real[1:] = random_state.normal(0, 1, len(omega) - 1)

    x_fft.imag[1:] = random_state.normal(0, 1, len(omega) - 1)

#    x_fft[1:] *=2*sigma*np.sqrt(tau/(1+omega[1:]**2*tau**2))
    x_fft[1:] *= np.sqrt(2)*sigma*np.sqrt(tau/(1+omega[1:]**2*tau**2))

    x_fft[1:] *= (1. / np.sqrt(2))


    if (not generate_complex) and (N % 2 == 0):
    
        x_fft.imag[-1] = 0

    if generate_complex:
        
        x = np.fft.ifft(x_fft)

    else:
        
        x = np.fft.irfft(x_fft, N)




    return x

