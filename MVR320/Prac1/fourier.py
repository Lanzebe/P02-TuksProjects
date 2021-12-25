import numpy as np
def fint(x,fs,IS=1,lc=0.0):

	""" fint(x,fs,IS=1,lc=0.0)

	fint      Integrates a time domain signal IS times in the frequency domain.

	Parameters
	----------
	x       Time signal to be integrated.  Can be 1D list or 1D Numpy array.
	fs      Sampling frequency [Hz]
	IS      Number of integration steps to be performed, default to 1.
	lc      Lower cut-off frequency [Hz], default to 0
	Returns
	-------
	y	Integrated signal
	t	Time array corresponding to the elements in y

	The frequency components of time signal x lying between 0 Hz and lc 
	is zeroed in the frequency domain.
	Furthermore 0-padding of the signal may be done to get a signal length of 2^n.
	This will greatly reduce the computational time needed for the FFT,
	and IFFT operations.
	Dynamic Systems Group/Centre for Asset Integrity Management
	HJ van Niekerk 1993-08-03 
	Modified by PS Heyns 1999-05-24; 2006-05-29
	Modified by DH Diamond 2014-08-22 """

	h = np.fft.rfft(np.array(x)) # Compute the fft of a real input	
	N = len(h)
	df = fs/float(len(x))
	freq = np.arange(N) * df # Get the frequencies
	cutoff_index = np.where(freq < lc)[0]
	if len(cutoff_index) != 0:
		cutoff_index = cutoff_index[-1]
		h[:cutoff_index+1] = np.zeros(len(h[:cutoff_index+1])) # Get the cutoff frequency
	omega = freq * 2 * np.pi
	omega[0] = 1
	g = h / (1j * omega)**IS	
	y = np.fft.irfft(g)
	return y, np.arange(len(y)) * 1/fs
