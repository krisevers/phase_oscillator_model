import numpy as np
from scipy.signal import savgol_filter

def noise(number_of_oscillators, simulation_time, dt):
    # %%        Noise        ####
    # Noise is important to be able to distinguish true synchrony from false synchrony
    # (= insufficient dephasing)
    noiseterm = np.zeros((number_of_oscillators,int(simulation_time/dt)))
    for ind in range(number_of_oscillators):
        # from Little et al. (2007), "Exploiting nonlinear recurrence and fractal scaling properties for voice
        # disorder detection", Biomed Eng Online, 6:23
        N       = (simulation_time/dt)                          # number of time steps
        alpha   = 1                                             
        N2      = np.floor(N/2)-1                                
        f       = np.transpose((np.arange(1,N2+1,1)))       
        A1      = 1/((N2+2)**alpha)                             
        A2      = 1/(f**(alpha/2))                              
        p2      = (np.random.random([int(N2),1])-0.5)*2*np.pi
        VW      = np.exp(1j*p2)[:,0]
        d2      = A2 * np.real(VW)
        d       = np.append(1, d2)
        d       = np.append(d, A1 )
        d       = np.append(d, np.flipud(np.conj(d2)) )
        d       = np.transpose(d)
        x       = np.real(np.fft.ifft(d))

        noiseterm[ind,:] = ((x-np.mean(x))/np.std(x))/50 # phase detuning of std = 0.02
    noiseterm = savgol_filter(noiseterm, 3, 2)

    return noiseterm
