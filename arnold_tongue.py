def arnold_tongue(C_vec, orientations, len_of_network, trials=20, show:bool=True):
    num_oscs = len_of_network ** 2
    at = np.zeros((len(orientations), len(C_vec)))
    for i in tqdm(range(len(C_vec))):
        for it in range(len(orientations)):
            bins = 0
            for k in range(trials):
                A, W, image, box, x_osc, y_osc, largest_contrasts, osc_of_interest, number_of_oscillators = stimulus(
                    number_of_oscillators=num_oscs, size=len_of_network, set='random', set_parameter=orientations[it], input_only=True,
                    bg=True, bg_fixed=True, bg_density=1, show=False)

                K, dist, alignment = connectivity(x_osc, y_osc, A, W, s_dist, s_align, spec_align, C_vec[i])

                initial_phase = np.random.rand(number_of_oscillators,20)*1*(2*np.pi) #Initial phases
                phases = np.zeros((number_of_oscillators,int(simulation_time/dt)))
                phases[:,0:20] = initial_phase

                noiseterm = noise(number_of_oscillators, simulation_time, dt)

                allcoh, allph = kuramoto(K, W, phases, noiseterm, simulation_time, dt, number_of_oscillators)

                path_oscs = np.where(W == np.max(W))
                for index, x in np.ndenumerate(alignment[path_oscs,path_oscs]):
                    bins += allcoh[index]

                bins = bins / trials
                at[it,i] = bins

    if show:
        plt.figure(figsize=(10,10))
        plt.imshow(at.T)
        plt.title('Arnold Tongue | s_dist = {}, s_align = {}, spec_align = {}'.format(s_dist, s_align, spec_align))
        plt.xlabel('orientatation')
        plt.ylabel('C')
        plt.xticks(np.arange(len(orientations)), orientations)
        plt.yticks(np.arange(len(C_vec)), np.around(C_vec, 4))
        plt.show()

    return at
