#!env python3


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.fftpack
import math
import subprocess
import wavio
import sys
# import sounddevice as sd
from optparse import OptionParser

import signal_env
from room_response_estimator import *

#######################################################################################################################

def Spectrum(s):
    Ftest = scipy.fftpack.fft( s )
    n = round(s.shape[0]/2)
    xf = np.linspace(0.0, 44100/2.0, n)
    return xf, 20*np.log10(np.abs(Ftest[0:n])) 

#######################################################################################################################
if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-m", "--measure-only", action="store_true",
                      dest="measure_only", default=False,
                      help="Measure without computing the deconvolution.")
    parser.add_option( "-r", "--reuse", action="store", 
                        type="string", dest="reuse_wav",
                        help="Use wav file with previous record instead actual playing and recording.")
    parser.add_option( "-d", "--duration", action="store", 
                        type="float", dest="duration",
                        default=10,
                        help="Duration of probe impulse.")
    parser.add_option( "-b", "--low-freq", action="store", 
                        type="float", dest="lowfreq",
                        default=100,
                        help="The lowest detected frequency [Hz].")
    parser.add_option( "-e", "--high-freq", action="store", 
                        type="float", dest="highfreq",
                        default=15000,
                        help="The highest frequency in probe impulse [Hz].")
    (options, args) = parser.parse_args()


    estimator = RoomResponseEstimator(options.duration, options.lowfreq, options.highfreq)

    if options.reuse_wav:
        ir_file = options.reuse_wav
    else:
        # Store probe signal in wav file.
        x = np.append(np.zeros(44100), estimator.probe_pulse)
        wavio.write("test_sweep.wav", x, 44100, sampwidth=2)

        # Play it and record microphones input simultaneously.
        reccommand = \
            "arecord -D plughw:1,0 -d {0} -f S16_LE -c2 -r44100 record_sweep.wav".format(
                int(options.duration + 5)).split(" ")
            # "rec -q --clobber -r 44100 -b 16 -D -c 2 record_sweep.wav trim 0 10".split(" ")
        prec = subprocess.Popen(reccommand)
        print(reccommand)

        playcommand = \
            "mplayer -volume 1 -really-quiet test_sweep.wav".split(" ")
        pplay = subprocess.Popen(playcommand)
        pplay.wait()
        prec.wait()

        ir_file = "record_sweep.wav"

    if options.measure_only:
        sys.exit(0)

    # Get the result of measurement from wav file.
    ir_fl = wavio.read( ir_file )
    # scale between -1 and 1, remove 2nd channel
    ir = ir_fl.data[0:,0]/math.pow(2.0, ir_fl.sampwidth*8-1)

    # Restore Room Response from the raw signal.
    room_response = estimator.estimate(ir)

    # Derive inverted room response for active room correction.
    Hi = fft(room_response)
    lmbd = 1e-2 # power spectral density of noise
    # Perform Weiner deconvolution.
    inv_room_response = np.real(ifft(np.conj(Hi)/(Hi*np.conj(Hi) + lmbd**2)))
    inv_room_response /= np.max(np.abs(inv_room_response))

    deconvolved_ir = fftconvolve(room_response, inv_room_response)
    deconvolved_sweep = fftconvolve(ir, inv_room_response)
    # save the deconvolved sweep
    deconvolved_wav = deconvolved_sweep * math.pow(2.0, ir_fl.sampwidth * 8 - 1)
    deconvolved_wav_16 = deconvolved_wav.astype(np.int16)
    wavio.write("deconvolved_sweep.wav", deconvolved_wav_16, 44100, sampwidth=2)
#######################################################################################################################
    plt.rcParams["figure.figsize"] = (6.4 * 2, 4.8 * 2)
    plt.subplot(321)
    plt.plot(room_response)
    plt.legend(["Measured Room Response"])
    x0 = np.argmax(room_response)
    plt.xlim([x0, x0+1024])
    plt.grid()

    plt.subplot(322)
    plt.plot(deconvolved_ir)
    plt.legend(["Deconvolved Room Response"])
    x0 = np.argmax(deconvolved_ir)
    plt.xlim([x0-100, x0+100])
    ymin, ymax = plt.ylim()
    plt.ylim( [ymin, ymax*1.6])
    plt.grid()

    plt.subplot(323)
    plt.plot(ir)
    plt.legend(["Recorded log-sine sweep"])
    plt.grid()

    ax = plt.subplot(324)
    plt.plot(*Spectrum(ir))
    plt.legend(["Spectrum of the record"])
    ax.set_xscale('log')
    ymin, ymax = plt.ylim()
    plt.ylim( [ymin, ymax*1.6])
    plt.xlim([10, 2e4])
    plt.grid()

    plt.subplot(325)
    plt.plot(deconvolved_sweep)
    plt.legend(["Deconvolved sine sweep"])
    plt.grid()

    ax = plt.subplot(326)
    plt.plot( *Spectrum(deconvolved_sweep))
    plt.legend(["Spectrum of deconvolved sine sweep" ])
    ax.set_xscale('log')
    ymin, ymax = plt.ylim()
    plt.ylim( [ymin, ymax*1.6])
    plt.xlim([10, 2e4])
    plt.grid()

    # plt.show()
    plt.savefig('measure_plt.png')
    # plt.close(fig)