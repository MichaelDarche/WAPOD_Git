#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 10:07:29 2025

@author: michael
"""

import wave
import numpy as np
import sys
import matplotlib.pyplot as plt

def wave2np(file):
    name=file+".wav"
    wav=wave.open(name,"rb")
    sampleRate=wav.getframerate()
    nsamples=wav.getnframes()
    raw=wav.readframes(-1)
    sampWidth=wav.getsampwidth()
    wav.close()
    taudio=nsamples/sampleRate
    raw=np.frombuffer(raw,dtype=np.int16)
    Time=np.linspace(0,taudio,num=nsamples)

    # np.savetxt("SEECUPP.txt", raw)
    return Time,raw,sampleRate,sampWidth,nsamples

def interpol(t,raw,dt):
    nt=int(t[-1]/dt)
    dt0=t[1]
    res=np.zeros(nt)
    for k in range(nt):
        ti=k*dt
        j=int(ti/dt0)
        ratio = (ti-dt0*j)/dt0
        res[k]=(1-ratio)*raw[j]+ratio*raw[j+1]
    return res
    
def np2wav(file,npdata,sR,sW,n):
    name=file+".wav"
    wav=wave.open(name,"wb")
    wav.setframerate(sR)
    wav.setnframes(n)
    wav.setnchannels(1)
    wav.setsampwidth(sW)
    data=npdata.tobytes()
    wav.writeframes(data)
    wav.close()
    # taudio=n/sR
    # Time=np.linspace(0,taudio,num=n)

    # np.savetxt("SEECUPP.txt", raw)
    return data

    
    
    
# sonwav=wave.open("Coq.wav","rb")
# sampleRate=sonwav.getframerate()
# nsamples=sonwav.getnframes()
# sW=sonwav.getsampwidth()
# raw=sonwav.readframes(-1)
# sonwav.close()

# taudio=nsamples/sampleRate

# raw2=np.frombuffer(raw,dtype=np.int16)


# Time=np.linspace(0,taudio,num=nsamples)


# sonwavmod=wave.open("CoqRev.wav","wb")
# sonwavmod.setnchannels(1)
# sonwavmod.setframerate(sampleRate)
# sonwavmod.setnframes(nsamples)
# sonwavmod.setsampwidth(sW)
# raw2=np.flip(raw2)
# data=raw2.tobytes()
# sonwavmod.writeframes(data)
# sonwavmod.close()
# np.savetxt("SEECUPP.txt", raw)