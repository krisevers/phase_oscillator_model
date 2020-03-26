def kuramoto(K, W, phases, noiseterm, simulation_time, dt, number_of_oscillators, order_parameter:bool=False):
    '''
    Kuramoto simulation. Simulates synchrony between oscillators.
    Takes the weights matrix (K) and the intrinsic frequency (W) as input.
    '''
    #np.disp('start simulation')
    for time in range(20,int(simulation_time/dt)):
        #if np.mod(time,500) == 0:
            #np.disp((str(time/(0.001/dt)) + 'ms of ' + str(simulation_time*1000) + 'ms'))
        for it in range(number_of_oscillators):
            interact        = (np.sin((phases[it,time-1] - phases[:,time-1])))
            phases[it,time] = phases[it,time-1] + (dt*W[it] + np.nansum(dt*K[:,it]*-interact)) + noiseterm[it,time]
    #np.disp('done simulation')

    ph      = np.transpose(np.mod(np.transpose(phases),2*np.pi))   # phases
    xx      = np.exp(1j*ph[:,99:])
    xxt     = np.transpose(xx)
    ab_real = np.real(xx).dot(np.real(xxt))
    ab_imag = np.imag(xx).dot(np.imag(xxt))
    ab      = ab_real + ab_imag*1J
    allcoh  = np.abs(ab)
    allcoh  = allcoh / np.max(allcoh)    # phase locking matrix
    allph   = np.angle(ab)               # phase relation matrix

    if order_parameter == True:
        itercoh = np.zeros((len(K),len(K),len(xxt)))
        iterph  = np.zeros((len(K),len(K),len(xxt)))
        allcoh  = np.zeros((len(K),len(K),len(xxt)))
        allph   = np.zeros((len(K),len(K),len(xxt)))
        for i in range(len(xxt)):
            ab_real = np.real(xx[:,i-1:i]).dot(np.real(xxt[i-1:i,:]))
            ab_imag = np.imag(xx[:,i-1:i]).dot(np.imag(xxt[i-1:i,:]))
            ab      = ab_real + ab_imag*1j
            itercoh[:,:,i]   = np.abs(ab)
            itercoh[:,:,i]   = itercoh[:,:,i] / np.max(itercoh[:,:,i])
            iterph[:,:,i]    = np.angle(ab)

            ab_real = np.real(xx[:,:i]).dot(np.real(xxt[:i,:]))
            ab_imag = np.imag(xx[:,:i]).dot(np.imag(xxt[:i,:]))
            ab      = ab_real + ab_imag*1j
            allcoh[:,:,i]   = np.abs(ab)
            allcoh[:,:,i]   = allcoh[:,:,i] / np.max(allcoh[:,:,i])
            allph[:,:,i]    = np.angle(ab)

        return allcoh, allph, itercoh, iterph

    else:       
        return allcoh, allph
