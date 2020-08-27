# ziaeeNN2020
-  This repository contains the codes for repuducing results and figures in the paper: 
-  Frequency-dependent organization of the brain’s functional network through delayed-interactions [J. Neural Networks](https://www.researchgate.net/publication/343780297_Frequency-dependent_organization_of_the_brain's_functional_network_through_delayed-interactions) 
- Abolfazl Ziaeemehr, Mina Zarei, Alireza Valizadeh, Claudio R. Mirasso

-  **Abstract**

*The structure of the brain network shows modularity at multiple spatial scales. The effect of the modular structure on the brain dynamics has been the focus of several studies in recent years but many aspects remain to be explored. For example, it is not well-known how the delays in the transmission of signals between the neurons and the brain regions, interact with the modular structure to determine the brain dynamics. In this paper, we show an important impact of the delays on the collective dynamics of the brain network with modular structure; that is, the degree of the synchrony between different brain regions is dependent on the frequency. In particular, we show that increasing the frequency the network transits from a global synchrony state to an asynchronous state, through a transition region over which the local synchrony inside the modules is stronger than the global synchrony. When the delays are dependent on the distance between the nodes, the modular structure of different spatial scales appears in the correlation matrix over different specific frequency bands, so that, finer spatial modular structure reveal in higher frequency bands. The results are justified by a simple theoretical argument and elaborated by simulations on several simplified modular networks and the connectome with different spatial resolutions.*

-  **How to use:**
```

$ cd codes/src
# set the parameters in run.py
$ make clean [or eradicate]
$ make 
$ python3 run.py
```

If you used the code please cite the paper:
- Ziaeemehr, A., Zarei, M., Valizadeh, A. and Mirasso, C.R., 2020. Frequency-dependent organization of the brain’s functional network through delayed-interactions. Neural Networks. https://doi.org/10.1016/j.neunet.2020.08.003
