'''
Small macro retrieving data from a gwf file a putting if into HR (Human-Readable) format

If injections are present they are also added to the output file


SV: 18/12/2023

Option:

->ifile     : input file name (gwf file)
->ofile     : output file name (ROOT file)
->fe        : sampling frequency of output data, in Hz (Default is 2048).
->d         : duration to extract, in s (Default is 100).


'''

import numpy as npy

try:
    import ROOT as root
except ImportError as e:
    print("ROOT is not installed, we will only produce a pickle file")
    pass  # module doesn't exist, deal with it.

import pickle
import os
import re
import math
from array import array
from virgotools import getChannel
from gwpy.timeseries import TimeSeries, TimeSeriesDict


def parse_cmd_line():
    import argparse
    """Parseur pour la commande gendata"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input","-i",help="Input file name",default=None)
    parser.add_argument("--output","-o",help="Input file name",default=None)
    parser.add_argument("-fe",help="Output signal frequency",type=float,default=2048)
    parser.add_argument("-d",help="Duration to extract (in s)",type=float,default=100)

    args = parser.parse_args()
    return args



'''
Main part of gendata.py
'''

def main():

    # Start by parsing the input options
    args = parse_cmd_line()

    filename=args.input
    print("Step 1: get the gwf input file info")

    # Quick'n'dirty file info
    parsename= "fileInfo.txt"
    os.system(f"FrDump -i {filename} -d 0  > {parsename}")

    filen = open(parsename, 'r')
    lines = filen.readlines()
    fileinfo=re.split(' |\t',lines[0])
    os.system(f"rm {parsename}")

    tinit    = float(fileinfo[1])
    duration = args.d
    fs       = args.fe

    print("Step 2: get the channel data and process it")
    # Quick'n'dirty strain info
    parsename= "injStuff.txt"
    os.system(f"FrDump -i {filename} -event -f {tinit} -l {duration} -d 3  > {parsename}")

    # Get the relevant channels
    # Hanford and Livingston strains, with and without injections

    testH_I = getChannel(filename,"H1:GDS-CALIB_STRAIN_CLEAN_AR_INJOFFLINE_4096Hz",tinit,tinit+duration) 
    testH   = getChannel(filename,"H1:GDS-CALIB_STRAIN_CLEAN_AR_4096Hz",tinit,tinit+duration) 
    testL_I = getChannel(filename,"L1:GDS-CALIB_STRAIN_CLEAN_AR_INJOFFLINE_4096Hz",tinit,tinit+duration)
    testL   = getChannel(filename,"L1:GDS-CALIB_STRAIN_CLEAN_AR_4096Hz",tinit,tinit+duration)

    # Get info as time serie
    TSH  = TimeSeries(data=testH.data, t0=testH.gps, sample_rate=testH.fsample,unit="", channel=testH.name)
    TSHI = TimeSeries(data=testH_I.data, t0=testH_I.gps, sample_rate=testH_I.fsample,unit="", channel=testH_I.name)
    TSL  = TimeSeries(data=testL.data, t0=testL.gps, sample_rate=testL.fsample,unit="", channel=testL.name)
    TSLI = TimeSeries(data=testL_I.data, t0=testL_I.gps, sample_rate=testL_I.fsample,unit="", channel=testL_I.name)

    # Resample properly
    TSH  = TSH.resample(fs)
    TSHI = TSHI.resample(fs)
    TSL  = TSL.resample(fs)
    TSLI = TSLI.resample(fs)

    delta_t=TSH.dt.value

    # Here we finally have the strain in array format
    Lh_of_t=TSL.data
    Lh_of_ti=TSLI.data
    Hh_of_t=TSH.data
    Hh_of_ti=TSHI.data
     
    # Spare some space 
    del testH_I,testH,testL_I,testL
    del TSH,TSHI,TSL,TSLI

    time     = array('d', [0])
    strainH  = array('d', [0])
    strainHI = array('d', [0])	
    strainL  = array('d', [0])
    strainLI = array('d', [0])

    lumiD    = array('d', [0])

    M1     = array('d', [0])
    M2     = array('d', [0])
    M1s    = array('d', [0])
    M2s    = array('d', [0])
    s1x    = array('d', [0])
    s1y    = array('d', [0])
    s1z    = array('d', [0])
    s2x    = array('d', [0])
    s2y    = array('d', [0])
    s2z    = array('d', [0])
    chieff = array('d', [0])
    chip   = array('d', [0])


    t_coal = array('d', [0])
    t_coalH = array('d', [0])
    t_coalL = array('d', [0])
    SNR_H  = array('d', [0])
    SNR_L  = array('d', [0])

    # Store everything into a ROOT file
    file  = root.TFile.Open(args.output, 'recreate')
    tree  = root.TTree("strain", "strain")
    tree3 = root.TTree("injections", "injections")

    tree.Branch("h_of_ti_H", strainHI,'strainHI/D')
    tree.Branch("h_of_t_H",   strainH,'strainH/D')
    tree.Branch("h_of_ti_L", strainLI,'strainLI/D')
    tree.Branch("h_of_t_L",   strainL,'strainL/D')
    tree.Branch("t",          time,  'time/D')

    tree3.Branch("mass1",     M1,     'M1/D')
    tree3.Branch("mass2",     M2,     'M2/D')
    tree3.Branch("mass1_s",   M1s,    'M1s/D')
    tree3.Branch("mass2_s",   M2s,    'M2s/D')
    tree3.Branch("tcoalH",    t_coalH,'t_coalH/D')
    tree3.Branch("tcoalL",    t_coalL,'t_coalL/D')
    tree3.Branch("tcoal",     t_coal, 't_coal/D')
    tree3.Branch("SNR_H",     SNR_H,  'SNR_H/D')
    tree3.Branch("SNR_L",     SNR_L,  'SNR_L/D')
    tree3.Branch("s1x",       s1x,    's1x/D')
    tree3.Branch("s1y",       s1y,    's1y/D')
    tree3.Branch("s1z",       s1z,    's1z/D')
    tree3.Branch("s2x",       s2x,    's2x/D')
    tree3.Branch("s2y",       s2y,    's2y/D')
    tree3.Branch("s2z",       s2z,    's2z/D')
    tree3.Branch("chi_eff",   chieff, 'chieff/D')
    tree3.Branch("chi_p",     chip,   'chip/D')
    tree3.Branch("Deff",      lumiD,  'lumiD/D')

    compt=0

    for h in Lh_of_t:
        if compt*delta_t>=duration:
            break
        time[0]=tinit+compt*delta_t
        
        # Those little checks are there to make sure we got gated data too
        if math.isnan(Lh_of_t[compt]):
            strainL[0]=0.
        else:
            strainL[0]=Lh_of_t[compt]

        if math.isnan(Lh_of_ti[compt]):
            strainLI[0]=0.
        else:
            strainLI[0]=Lh_of_ti[compt]

        if math.isnan(Hh_of_t[compt]):
            strainH[0]=0.
        else:
            strainH[0]=Hh_of_t[compt]

        if math.isnan(Hh_of_ti[compt]):
            strainHI[0]=0.
        else:
            strainHI[0]=Hh_of_ti[compt]

        tree.Fill()
        compt+=1


    print("Step 3: get the injection data")

    # The dirty little parser
    filen = open(parsename, 'r')
    lines = filen.readlines()
    injection=[]
    for line in lines:
        if "SimEvent:" in line:
            if "FrSimEvent:" in line:
                continue    
            params=[]
            words=re.split(' |\t',line)
            for word in words:
                if "=" in word:
                     vals=word.split("=")
                     params.append((vals[0],float(vals[1])))
            injection.append(params)
        else: 
            continue
    
    filen.close()

#
# Params defined here:
# https://git.ligo.org/RatesAndPopulations/o4-plan-investigations/-/blob/main/o4-injection-proposal/o4-file-formats.md#deterministic-derived-or-fixed-properties-of-the-injection
#


    for inj in injection:
        for pars in inj:
            if pars[0]=='time_H':
                t_coalH[0]=pars[1]
            if pars[0]=='time_L':
                t_coalL[0]=pars[1]
            if pars[0]=='time_geocenter':
                t_coal[0]=pars[1]
            if pars[0]=='mass1_detector':
                M1[0]=pars[1]
            if pars[0]=='mass2_detector':
                M2[0]=pars[1]
            if pars[0]=='mass1_source':
                M1s[0]=pars[1]
            if pars[0]=='mass2_source':
                M2s[0]=pars[1]
            if pars[0]=='observed_snr_L':
                SNR_L[0]=pars[1]
            if pars[0]=='observed_snr_H':
                SNR_H[0]=pars[1]
            if pars[0]=='spin1x':
                s1x[0]=pars[1]
            if pars[0]=='spin1y':
                s1y[0]=pars[1]
            if pars[0]=='spin1z':
                s1z[0]=pars[1]
            if pars[0]=='spin2x':
                s2x[0]=pars[1]
            if pars[0]=='spin2y':
                s2y[0]=pars[1]
            if pars[0]=='spin2z':
                s2z[0]=pars[1]
            if pars[0]=='chi_eff':
                chieff[0]=pars[1]
            if pars[0]=='chi_p':
                chip[0]=pars[1]
            if pars[0]=='luminosity_distance':
                lumiD[0]=pars[1]
        tree3.Fill()

    file.Write()
    file.Close()

    # Clean the dirty stuff
    os.system(f"rm {parsename}")


############################################################################################################################################################################################
if __name__ == "__main__":
    main()
