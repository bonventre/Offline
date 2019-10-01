import numpy as np
import scipy.signal as sigpy
from math import *

def omega2freq(omega):
  return -1*omega/(2*np.pi)

def freq2omega(freq):
  return -1*freq*2*np.pi

def digital_coefficients(poles,zeros,fs):
  k = 1.
  poles = map(freq2omega,poles)
  zeros = map(freq2omega,zeros)
  analog_coefficients = sigpy.zpk2tf(zeros,poles,k)
  return sigpy.bilinear(analog_coefficients[0],analog_coefficients[1],fs)

def frequency_response(poles,zeros,fs):
  dc = digital_coefficients(poles,zeros,fs)
  w,h = sigpy.freqz(dc[0],dc[1],2**20)
  w = w * fs / (2*np.pi)
  return w,h

def frequency_response_analog(poles,zeros):
  k = 1.
  poles = map(freq2omega,poles)
  zeros = map(freq2omega,zeros)
  ac = sigpy.zpk2tf(zeros,poles,k)
  w,h = sigpy.freqs(ac[0],ac[1],2**20)
  w = w / (2*np.pi)
  return w,h



def apply_filter(poles,zeros,fs,input_pulse):
  dc = digital_coefficients(poles,zeros,fs)
  return sigpy.lfilter(dc[0],dc[1],input_pulse)

def apply_filter_from_dc(dc,input_pulse):
  return sigpy.lfilter(dc[0],dc[1],input_pulse)
