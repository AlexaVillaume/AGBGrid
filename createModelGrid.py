'''
Description: Generate the grid of models for the dusty envelopes around TP-AGB stars

Author: Alexa Villaume

Date: 10/8/13 - 10/14/13
'''

import sys, os, commands, math, string
import glob
import numpy as np
from subprocess import call

#def createOSpec(input):

def createCSpec(input):
	spectrum = np.loadtxt(input)
	wave = spectrum[:,0]
	flux = spectrum[:,2]

	# Convert the Aringer wavelength grid to microns
	lam_mu = []
	for i in range(len(wave)):
		lam_mu.append(wave[i]*0.0001)

	out = np.column_stack((lam_mu, flux))
	return out

def generateInput(tau, co, spectra):
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
		dusttemp = '1063.35'
		SC = '0.90'
		SiC_Pg = '0.10'
	else:
		dusttemp = '1331.47'
		Sil_Ow = '1.00'


	f = open("temp.inp", 'w')
	f.write("""
  I PHYSICAL PARAMETERS
     1) External radiation:
        Spectrum = 4
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

	f.close()

def makeCgrid(i, j):
	call(['./dusty'])
	call(['mv', 'temp.stb', 'CGrid/temp_'+'teff'+i[9:13]+'_'+'tau'+str(j)])
	return np.loadtxt('CGrid/temp_'+'teff'+i[9:13]+'_'+'tau'+str(j))

def makeOgrid(i, j):
	call(['./dusty'])
	call(['mv', 'temp.stb', 'OGrid/temp_'+'teff'+i[9:13]+'_'+'tau'+str(j)])

# Parse the DUSTY output into a single file
def PutItAllTogether(file, Teff, Tau, model, count, count2):
	if count == 0 and count2 == 0:
		wave = model[:,0]
		file.write('# Wavelength Grid:' + ' ' + str(len(wave)) + '\n')
		for i in range(len(wave)):
			file.write('%10.2e' % wave[i]*1e4)
		file.write('\n')
	flux = []
	fTot = model[:,1]
	fInp = model[:,5]
	file.write('#' + ' ' +Teff + '%10.2f' % math.log(Tau, 10) + '\n')
	for j in range(len(fTot)):
        tmp = fTot[j]/fInp[j]
        if tmp == float('NaN') or tmp == float('Inf'):
            flux.append(0)
        else:
		    flux.append(fTot[j]/fInp[j])
		file.write('%10.2e' % flux[j])
	file.write('\n')

def main():
    tau_min = 1e-3
    tau_max = 50.0
    tau = np.linspace(tau_min, tau_max, num=50, endpoint=True)
    taulog = []
    i = 1
    while i <= len(tau):
        q2 = math.exp(math.log(tau_max/tau_min)/(len(tau)-1))
        taulog.append(tau_min*q2**(i-1.0))
        i+=1

    co = 0
    while co <= 1:
        if co == 0:
	        f = open("OGrid.txt", "w")
            spec_in = (glob.glob('mxcom*'))
	        for i in range(len(spec_in)):
		        for j in range(len(taulog)):
			        out = createOSpec(spec_in[i])
			        np.savetxt("dustyin", out)
			        generateInput(str(tau[j]), co, "dustyin")
			        model = makeOgrid(spec_in[i], taulog[j])
			        PutItAllTogether(f, spec_in[i][9:13], taulog[j], model, i, j)
        else:
	        f = open("CGrid.txt", "w")
            spec_in = (glob.glob('mxcom*'))
	        for i in range(len(spec_in)):
		        for j in range(len(taulog)):
			        out = createCSpec(spec_in[i])
			        np.savetxt("dustyin", out)
			        generateInput(str(tau[j]), co, "dustyin")
				    model = makeCgrid(spec_in[i], taulog[j])
			        PutItAllTogether(f, spec_in[i][9:13], taulog[j], model, i, j)

if __name__ == "__main__":
	main()
