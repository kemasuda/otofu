import numpy as np

#%% jasmine noise model in ppm as a function of Hw=0.7J+0.3H for 5min
def noise_model(hw, hwlim=[7.756708395177219, 12.756708395177219], coeff=[2.73659218e+00, -9.95974284e+01,  1.37688923e+03, -8.47717270e+03,1.95861269e+04]):
    m = np.poly1d(coeff)
    if np.sum((hw<hwlim[0])|(hwlim[1]<hw)):
        print ('Hw outside the valid range:', hwlim)
    return m(hw)

#%% orbital period corresponding to the classical HZ around M dwarfs
def HZperiod(mass, s=1, z=np.poly1d([ 0.40127856,  1.54449068,  4.0068301 , -1.17924394])):
    return 12*s**(-3./4.)/np.sqrt(mass/0.15)*(np.exp(z(np.log(mass)))/3e-3)**(3./4.)
