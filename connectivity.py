def connectivity(x_osc, y_osc, A, W, s_dist, C):
    '''
    input should be a matrix with the dimensions of the 2D network and filled with zeros except for the locations
    where there is an signal with a certain angle.
    '''

    # distance weights
    dist        = np.zeros((len(x_osc),len(x_osc)))
    for i in range(len(x_osc)):
        dist[i]        = np.sqrt((x_osc[i]-x_osc[:])**2+(y_osc[i]-y_osc[:])**2)
        dist[i]        = np.maximum(0,dist[i])
        dist[i,i]      = 0    # avoid connections with itself

    # connectivity matrix
    K = C * np.exp(-dist/s_dist)
    np.fill_diagonal(K, 0)
    return K
