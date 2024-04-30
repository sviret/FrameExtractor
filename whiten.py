'''
Small macro whitening real data

The data is retrieved from a ROOT file produced with the macro gwf_extractor.py
and is whitened using PSDs computed on data

Real data can pass to 0 in some case (detector issue,...). This sharp discontinuity can induce big instabilities in the PSD
calculation, therefore the signal at the transition is smooth with a Blackman window. The size of the window can be tuned.



SV: 18/12/2023

Options:

->ifile     : input file name
->ofile     : output file name
->fe        : sampling frequency of output data (Default is 2048).

'''

import numpy as npy
import scipy
import pickle
import csv
import os
import sys
import math
import time
import ROOT as root
import matplotlib.pyplot as plt
from array import array
from scipy.stats import norm
from scipy import signal
from scipy import ndimage


def parse_cmd_line():
    import argparse
    """Parseur pour la commande gendata"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input","-i",help="Input file name injections",default=None)
    parser.add_argument("--output","-o",help="Input file name",default=None)
    parser.add_argument("-fe",help="Output signal frequency",type=float,default=2048)
    parser.add_argument("-tb",help="Smoothing window duration (in s)",type=float,default=1)
    parser.add_argument("--length","-l",help="PSD total duration (in s)",type=float,default=50)
    parser.add_argument("-fmin",help="Fréquence minimale visible par le détecteur",type=float,default=20)
    parser.add_argument("-fmax",help="Fréquence maximale visible par le détecteur",type=float,default=1000)
    
    args = parser.parse_args()
    return args



'''
Main part of gendata.py
'''

def main():

    # Start by parsing the input options

    args = parse_cmd_line()

    __inputfile=args.input
    __outputfile=args.output
    
    __fe=args.fe
    __fmin=args.fmin
    __fmax=args.fmax
    __Tblack=args.tb
    __PSD_lat=args.length
  

    #Step 1 retrieve the info processed with gwf_extractor.py
    
    fileD  = root.TFile.Open(__inputfile,"READ")
    treeD  = fileD.Get("strain")
    treeI  = fileD.Get("injections")


    _brutePSDL=[]
    _brutePSDH=[]
    compt=0
    PSD_length=int(__PSD_lat*__fe)

    h_of_t_H=[]
    h_of_ti_H=[]
    h_of_t_L=[]
    h_of_ti_L=[]
    comptH=0
    comptL=0
    PSDsH=[]
    PSDsL=[]

    idx=0
    nentries=treeD.GetEntries()  # Number of samples
    ninj=treeI.GetEntries()  # Number of injections


    # From there one can compte the number of PSDs which will be processed
    
    __N=PSD_length
    __delta_t=1/__fe        # Time step
    __delta_f=__fe/__N      # Frequency step
    __norm=npy.sqrt(__N)
    x=npy.zeros(__N, dtype=npy.float64)
    y=npy.zeros(__N, dtype=npy.float64)
    ifmax=int(min(__fmax,__fe/2)/__delta_f)
    ifmin=int(__fmin/__delta_f)

    one_to_zero_H=[]
    one_to_zero_L=[]
    zero_to_one_H=[]
    zero_to_one_L=[]

    prevH=-1
    prevL=-1

    # Loop over the strain data to compude the PSDs along
    # time. We use the data w/o injections here of course
    
    for entry in treeD:

        valH=entry.h_of_t_H  # Hanford
        valL=entry.h_of_t_L  # Livingston

        # Look for possible transition

        # 0->1, start a new data stream
        if valH!=0:
            x[comptH]=valH
            comptH+=1
            if idx>0 and prevH==0:
                zero_to_one_H.append(idx)
        
        if valL!=0:
            y[comptL]=valL
            comptL+=1
            if idx>0 and prevL==0:
                zero_to_one_L.append(idx)

        # If 1->0 then reset the data vector
        if valH==0 and idx>0 and prevH!=0:
            comptH=0
            x[:]=0.
            one_to_zero_H.append(idx)
            
        if valL==0 and idx>0 and prevL!=0:
            comptL=0
            y[:]=0.
            one_to_zero_L.append(idx)

        # If enough data was collected, one compute the
        # PSD, using the Welch method
        # And reset the stream

        # Hanford
        if comptH==PSD_length-1:
            _,PSD_H=signal.welch(x,fs=__fe,nperseg=__fe)
            comptH=0
            x[:]=0.
            
            factor=float((__N//2+1)/len(PSD_H))
            __PSDH=npy.ones(__N, dtype=npy.float64)
            __PSDH[ifmin:ifmax]=npy.abs(ndimage.zoom(PSD_H,factor))[ifmin:ifmax]
            __PSDH[-1:-__N//2:-1]=__PSDH[1:__N//2]
            
            PSDsH.append((idx,__PSDH,PSD_H))

        # Livingston
        if comptL==PSD_length-1:
            _,PSD_L=signal.welch(y,fs=__fe,nperseg=__fe)
            comptL=0
            y[:]=0.
        
            factor=float((__N//2+1)/len(PSD_L))
            __PSDL=npy.ones(__N, dtype=npy.float64)
            __PSDL[ifmin:ifmax]=npy.abs(ndimage.zoom(PSD_L,factor))[ifmin:ifmax]
            __PSDL[-1:-__N//2:-1]=__PSDL[1:__N//2]
        
            PSDsL.append((idx,__PSDL,PSD_L))

        
        idx+=1
        prevH=valH
        prevL=valL
        if idx%100000==0:
            print(idx,nentries)

    #
    # Loop over data is now finished, PSDs have been computed
    #

    step=int(0.5*PSD_length)
    nblocks=int(nentries/step)

    print("0 to 1 transition in Livingston data:",zero_to_one_L)
    print("0 to 1 transition in Hanford data:",zero_to_one_H)
    print("1 to 0 transition in Livingston data:",one_to_zero_L)
    print("1 to 0 transition in Hanford data:",one_to_zero_H)

    print("Need to smooth those transitions before the whitening")
    print("Use Blackman windows for that")

    h_of_t_wH   = array('d', [0])
    h_of_ti_wH  = array('d', [0])
    h_of_t_wL   = array('d', [0])
    h_of_ti_wL  = array('d', [0])
    t           = array('d', [0])
    
    psdh        = array('d', [0])
    psdl        = array('d', [0])
    rank        = array('d', [0])
    
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

    file2 = root.TFile.Open(__outputfile, 'recreate')
    treep  = root.TTree("white", "white")
    treeh = root.TTree("psds_H", "psds_H")
    treel = root.TTree("psds_L", "psds_L")
    tree3 = root.TTree("injections", "injections")
       
    treep.Branch("h_of_t_white_H",   h_of_t_wH, 'h_of_t_wH/D')
    treep.Branch("h_of_ti_white_H",  h_of_ti_wH, 'h_of_ti_wH/D')
    treep.Branch("h_of_t_white_L",   h_of_t_wL, 'h_of_t_wL/D')
    treep.Branch("h_of_ti_white_L",  h_of_ti_wL, 'h_of_ti_wL/D')
    treep.Branch("time",             t,'t/D')

    treeh.Branch("PSD_H",   psdh, 'psdh/D')
    treel.Branch("PSD_L",   psdl, 'psdl/D')
    treeh.Branch("rank",   rank, 'rank/D')
    treel.Branch("rank",   rank, 'rank/D')

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

    
    PSDL=npy.ones(__N, dtype=npy.float64)
    PSDH=npy.ones(__N, dtype=npy.float64)

    PSDL=PSDsL[0][1]
    PSDH=PSDsH[0][1]

    idxH=[]
    idxL=[]

    # We store the PSD info in the output file

    rk=0
    for psd in PSDsL:
        idxL.append(psd[0])
        field=psd[2]
        for val in field:
            psdl[0]=val
            rank[0]=rk
            treel.Fill()
        rk+=1
    
    rk=0
    for psd in PSDsH:
        idxH.append(psd[0])
        field=psd[2]
        for val in field:
            psdh[0]=val
            rank[0]=rk
            treeh.Fill()
        rk+=1

    # Finally time for the whitening
    # Pick up the data for every block and do the whitening of the block

    for i in range(nblocks):

        if (i*step+PSD_length)>nentries: # End of file
            break

        start_i=i*step
        stop_i=i*step+PSD_length
    
        SH=npy.zeros(__N, dtype=npy.float64)
        SiH=npy.zeros(__N, dtype=npy.float64)
        SL=npy.zeros(__N, dtype=npy.float64)
        SiL=npy.zeros(__N, dtype=npy.float64)

        # Get the data
        for index in range(start_i,stop_i):
            treeD.GetEntry(index)
            val=index-start_i
            SH[val]=entry.h_of_t_H
            SL[val]=entry.h_of_t_L

            SiH[val]=entry.h_of_ti_H
            SiL[val]=entry.h_of_ti_L
            
    
        # Smooth the internal discontinuities (if applicable)
        for risingedge in zero_to_one_L:
            idx=risingedge-start_i
            if 0<idx<__N:
                iwmax=int(min(idx+__fe/2,__N-1))
                iwmin=int(idx)
                SL[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[:iwmax-iwmin]
                SiL[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[:iwmax-iwmin]
   
        for fallingedge in one_to_zero_L:
            idx=fallingedge-start_i
            if 0<idx<__N:
                iwmax=int(idx)
                iwmin=int(max(0,idx-__fe/2))
                SL[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[iwmax-iwmin:]
                SiL[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[iwmax-iwmin:]

        for risingedge in zero_to_one_H:
            idx=risingedge-start_i
            if 0<idx<__N:
                iwmax=int(min(idx+__fe/2,__N-1))
                iwmin=int(idx)
                SH[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[:iwmax-iwmin]
                SiH[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[:iwmax-iwmin]
   
        for fallingedge in one_to_zero_H:
            idx=fallingedge-start_i
            if 0<idx<__N:
                iwmax=int(idx)
                iwmin=int(max(0,idx-__fe/2))
                SH[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[iwmax-iwmin:]
                SiH[iwmin:iwmax]*=npy.blackman((iwmax-iwmin)*2)[iwmax-iwmin:]
            
        # Use the PSD closest to the current block
        compt=0
        for idx in idxL:
            if idx<i*step:
                break
            compt+=1
        if compt!=0:
            PSDL=PSDsL[min(compt,len(idxL)-1)][1]
        compt=0
        for idx in idxH:
            if idx<i*step:
                break
            compt+=1
        if compt!=0:
            PSDH=PSDsH[min(compt,len(idxH)-1)][1]

    
        print("Dealing with block ",i,"/",nblocks)

        __h_of_t_wH = npy.zeros(__N, dtype=npy.float64)
        __h_of_ti_wH = npy.zeros(__N, dtype=npy.float64)
        __h_of_t_wL = npy.zeros(__N, dtype=npy.float64)
        __h_of_ti_wL = npy.zeros(__N, dtype=npy.float64)

        __SfL=npy.zeros(__N, dtype=complex)
        __SfiL=npy.zeros(__N, dtype=complex)
        __SfH=npy.zeros(__N, dtype=complex)
        __SfiH=npy.zeros(__N, dtype=complex)

        iwmax_i=int(__Tblack/__delta_t)
        iwmin_i=0

        iwmax_f=__N
        iwmin_f=__N-int(__Tblack/__delta_t)
    
        # Smooth the block edges
        SL[iwmin_i:iwmax_i]*=npy.blackman((iwmax_i-iwmin_i)*2)[:iwmax_i-iwmin_i]
        SL[iwmin_f:iwmax_f]*=npy.blackman((iwmax_f-iwmin_f)*2)[iwmax_f-iwmin_f:]
        SiL[iwmin_i:iwmax_i]*=npy.blackman((iwmax_i-iwmin_i)*2)[:iwmax_i-iwmin_i]
        SiL[iwmin_f:iwmax_f]*=npy.blackman((iwmax_f-iwmin_f)*2)[iwmax_f-iwmin_f:]
        SH[iwmin_i:iwmax_i]*=npy.blackman((iwmax_i-iwmin_i)*2)[:iwmax_i-iwmin_i]
        SH[iwmin_f:iwmax_f]*=npy.blackman((iwmax_f-iwmin_f)*2)[iwmax_f-iwmin_f:]
        SiH[iwmin_i:iwmax_i]*=npy.blackman((iwmax_i-iwmin_i)*2)[:iwmax_i-iwmin_i]
        SiH[iwmin_f:iwmax_f]*=npy.blackman((iwmax_f-iwmin_f)*2)[iwmax_f-iwmin_f:]

        # FFT of the blocks
        __SfL=npy.fft.fft(SL,norm='ortho')
        __SfiL=npy.fft.fft(SiL,norm='ortho')
        __SfH=npy.fft.fft(SH,norm='ortho')
        __SfiH=npy.fft.fft(SiH,norm='ortho')

        # Normalisation and iFFT to go back to time domain
        __h_of_t_wH = npy.fft.ifft(__SfH[:]/npy.sqrt(__PSDH*__N*__delta_f/2.),norm='ortho').real
        __h_of_ti_wH = npy.fft.ifft(__SfiH[:]/npy.sqrt(__PSDH*__N*__delta_f/2.),norm='ortho').real
        __h_of_t_wL = npy.fft.ifft(__SfL[:]/npy.sqrt(__PSDL*__N*__delta_f/2.),norm='ortho').real
        __h_of_ti_wL = npy.fft.ifft(__SfiL[:]/npy.sqrt(__PSDL*__N*__delta_f/2.),norm='ortho').real

        # Store the relevant info in the output data file
        for j in range(int(0.25*__N),int(0.75*__N)):

            h_of_t_wH[0]=__h_of_t_wH[j]
            h_of_ti_wH[0]=__h_of_ti_wH[j]
            h_of_t_wL[0]=__h_of_t_wL[j]
            h_of_ti_wL[0]=__h_of_ti_wL[j]
            t[0]=(j+i*step)*__delta_t
            treep.Fill()
    

        # Clean
        del __SfL,__SfiL,__SfH,__SfiH
        del __h_of_t_wL,__h_of_ti_wH,__h_of_t_wH,__h_of_ti_wL
        del SH,SL,SiH,SiL
    



    # Copy the injections
    for entry in treeI:

        t_coalH[0]=entry.tcoalH
        t_coalL[0]=entry.tcoalL
        t_coal[0]=entry.tcoal
        M1[0]=entry.mass1
        M2[0]=entry.mass2
        M1s[0]=entry.mass1_s
        M2s[0]=entry.mass2_s
        SNR_L[0]=entry.SNR_L
        SNR_H[0]=entry.SNR_H
        s1x[0]=entry.s1x
        s1y[0]=entry.s1y
        s1z[0]=entry.s1z
        s2x[0]=entry.s2x
        s2y[0]=entry.s2y
        s2z[0]=entry.s2z
        chieff[0]=entry.chi_eff
        chip[0]=entry.chi_p
        lumiD[0]=entry.Deff
        tree3.Fill()
    
    fileD.Close()
    file2.Write()
    file2.Close()


############################################################################################################################################################################################
if __name__ == "__main__":
    main()
