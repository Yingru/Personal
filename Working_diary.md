# Working_diary

[TOC]

## 2018-07-30

- thermal boltzmann distribution
```python
def Boltzmann_E(E, temp, M):
    dfdE = np.exp(-E/temp) * np.sqrt(E**2 - M**2) * E
    normed = dfdE.sum()*(E[1] - E[0])
    return dfdE/normed
    
def Boltzmann_p(p, temp, M):
    E = np.sqrt(p**2 + M**2)
    dNdp = p/E * dNdE
    normed = dfdp.sum() * (p[1] - p[0])
    return dNdp/normed

# find the effective temperature of a system with equilibirum distribution
# being the boltzmann distribution
def effective_temperature(p0_list):
    from scipy.optmize import curve_fit
    dNdp0, bins = np.histogram(p0_list, bins=50, normed=True)
    p0_bins = 0.5*(bins[1:] + bins[:-1])
    def func(E, temp):
        return thermal_dfdE(E, temp, 1.3)
        
    popt, pcov = curve_fit(func, p0_bins, dNdp0)
    return popt  # array of parameters
```