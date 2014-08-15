'''
Description: Generate the grid of models for the dusty envelopes around TP-AGB stars
'''

import sys
import os
import commands
import math
import string
import glob
from subprocess import call
import numpy as np


"""
Using Aringer input spectra for the C-rich grid and
Basel input spectra for the O-rich grid.

DUSTY requires a wavelength grid in microns as input.
"""
def create_o_spectrum(input):
    flux = np.loadtxt(input)
    wave = np.loadtxt("basel.lambda")
    lam_mu = map(lambda line: line*1e-4, wave)
    out = np.column_stack((lam_mu, flux))
    return out

def create_c_spectrum(input):
    spectrum = np.loadtxt(input)
    wave = spectrum[:,0]
    flux = spectrum[:,2]
    lam_mu = map(lambda line: line*1e-4, wave)
    out = np.column_stack((lam_mu, flux))
    return out

def generate_input(tau, co, spectra, cold_sic):
    grainsize = '0.1'

    Sil_Ow = '0.00'
    Sil_Oc = '0.00'
    Sil_DL = '0.00'
    grf_DL = '0.00'
    amC_Hn = '0.00'
    SiC_Pg = '0.00'
    SC 	   = '0.00'
    SSC	   = '0.00'
    SSW    = '0.00'
    if co >= 1:
        spectrum = '4'
        dusttemp = '1100'
        SC = '0.90'
        SiC_Pg = '0.10'
    else:
        spectrum = '6'
        dusttemp = '700'
        if tau < 3:
            Sil_Ow = '1.00'
            Sil_Oc = '0.00'
        else:
            Sil_Ow = str(1 - cold_sic)
            Sil_Oc = str(cold_sic)
    with open("temp.inp", "w") as f:
	    f.write("""
  I PHYSICAL PARAMETERS
     1) External radiation:
        Spectrum = """+spectrum+"""
        """+spectra+"""

     2) Dust Properties

       2.1 Chemical composition
           optical properties index = 2
           Abundances for supported grain types:
               Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg
           x = """+Sil_Ow+"""  """+Sil_Oc+"""  """+Sil_DL+""" """+grf_DL+""" """+amC_Hn+""" """+SiC_Pg+"""
           Number of additional components = 3, properties listed in files
                    Suh_carbon.nk
                    Suh_silicate_c.nk
                    Suh_silicate_w.nk
           Abundaces for these components = """+SC+""", """+SSC+""", """+SSW+"""
       2.2 Grain size distribution

        - size distribution = 2 % Single size grain distribution
        - q = 0, a(min) = """+grainsize+""" micron, a(max) = """+grainsize+""" micron

       2.3 Dust temperature on inner boundary:

        - temperature = """+dusttemp+""" K


     3)Exponentially decreasing density distribution
        density type = 1; N = 1; Y = 1000; p = 2


     4) Optical Depth
        - grid type = 1                    % linear grid
        - lambda0 = 1.0 micron            % optical depth specified
        - tau(min) = """+tau+"""; tau(max) = """+tau+"""  % for the visual wavelength
        - number of models = 1
  ----------------------------------------------------------------------

  II NUMERICS

     - accuracy for flux conservation = 0.05

  ----------------------------------------------------------------------
  III OUTPUT PARAMETERS

    The separate flag 'verbose' controlls printout of messages to screen.
    If set to 0, there will be no screen output; if set to 1 - only minimal
    messages are printed out; if set to 2 - there will be more detailed
    screen output, in case the user would like to trace execution problems.
    The flags governing file production are as follows:
    If flag.eq.0 the particular file(s) is not produced. If flag.eq.1
    all model results are in corresponding files with extensions 'spp'
    (sp.properties), 'stb' (spectra), 'itb' (images and visibilities,
    if chosen), 'rtb' (radial profiles) and 'mtb' (messages).  If
    flag.eq.2 each model result is in a separate corresponding file,
    with visibilities contained in fname.i##. If the images flag.eq.3
    the visibilities will be in separate files fname.v## (the flag for
    visibilities has to be the same as for images).


        FILE DESCRIPTION                               FLAG
       ------------------------------------------------------------
       - verbosity flag;                               verbose = 1
       - properties of emerging spectra;             fname.spp = 1
       - detailed spectra for each model;           fname.s### = 1
       - images at specified wavelengths;           fname.i### = 0
       - radial profiles for each model;            fname.r### = 0
       - detailed run-time messages;                fname.m### = 0
       -------------------------------------------------------------

  The end of the input parameters listing.""")

def make_c_grid(temp, tau):
	call(['./dusty.exe'])
	call(['mv', 'temp.stb', 'CGrid/temp_'+'teff'+temp[9:13]+'_'+'tau'+str(tau)])
	return np.loadtxt('CGrid/temp_'+'teff'+temp[9:13]+'_'+'tau'+str(tau))

def make_o_grid(temp, tau):
	call(['./dusty.exe'])
	call(['mv', 'temp.stb', 'OGrid/temp_'+'teff'+temp[18:22]+'_'+'tau'+str(tau)])
	return np.loadtxt('OGrid/temp_'+'teff'+temp[18:22]+'_'+'tau'+str(tau))

# Parse the DUSTY output into a single file
def put_wave(file, model):
    file.write( str(len(model[:,0])) + '\n')
    wave = map(lambda value: value*1e4, model[:,0])
    file.write(''.join(['%10.2e' % value for value in  wave]))
    file.write('\n')

def put_all_together(file, Teff, Tau, model):
	flux = []
	fTot = model[:,1]
	fInp = model[:,5]
	file.write(Teff + '%10.2f' % math.log(Tau, 10) + '\n')

	for j in range(len(fTot)):
		ratio = fTot[j]/fInp[j]
		flux.append(ratio)
		if ratio == float('inf') or math.isnan(ratio):
			flux[j] = 0
		file.write('%10.2e' % flux[j])
	file.write('\n')

def main():
    if os.path.isdir('CGrid') == False:
        os.mkdir('CGrid')
    if os.path.isdir('OGrid') == False:
        os.mkdir('OGrid')

    # Get log spaced tau values
    tau_min = 1e-3
    tau_max = 50.0
    tau = np.linspace(tau_min, tau_max, num=50, endpoint=True)
    TauFidLog = []
    i = 1
    while i <= len(tau):
        q2 = math.exp(math.log(tau_max/tau_min)/(len(tau)-1))
        TauFidLog.append(round(tau_min*q2**(i-1.0), 5))
        i += 1

    """
    Cold and warm silicate switching for the O-rich grid.
    For higher values of tau want more cold silicate dust.
    """
    min_value = math.atan(tau_min)
    max_value = math.atan(tau_max)
    # The percent of cold silicates per tau value
    cold_sic = map(lambda value: 1 - ((max_value - math.atan(value))/max_value), TauFidLog)

    for i in range(len(TauFidLog)):
            print TauFidLog[i], cold_sic[i]

    whichgrid = 1
    while whichgrid <= 1:
        if whichgrid == 0:
            with open("CGrid.txt", "w") as gridfile:
                spec_in = glob.glob('mxcom*')
                for i, spectrum in enumerate(spec_in):
                    for j, tau in enumerate(TauFidLog):
                        print "Tau is", tau
                        np.savetxt("dustyin", create_c_spectrum(spectrum))
                        generate_input(str(tau), 1.1, "dustyin", cold_sic[j])
                        model = make_c_grid(spectrum, tau)
                        if i == 0 and j == 0:
                            put_wave(gridfile, model)
                        put_all_together(gridfile, spectrum[9:13], tau, model)
        if whichgrid == 1:
            with open("OGrid.txt", "w") as gridfile:
                spec_in = (glob.glob('z0.0190_logg0_temp*'))
                for i, spectrum in enumerate(spec_in):
                    for j, tau in enumerate(TauFidLog):
                        print "Tau is", tau
                        np.savetxt("dustyin", create_o_spectrum(spectrum))
                        generate_input(str(tau), 0.9, "dustyin", cold_sic[j])
                        model = make_o_grid(spectrum, tau)
                        if i == 0 and j == 0:
                            put_wave(gridfile, model)
                        put_all_together(gridfile, spectrum[18:], tau, model)
	whichgrid+=1

if __name__ == "__main__":
	main()
