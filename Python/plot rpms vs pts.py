import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

fs = 48000 # sampling frequency
f = 12000 # source frequency
fr = 1000 # rotation frequency

# sine panning curve
sine_gain = np.sin(np.linspace(1.5*np.pi, 3.5*np.pi, 512))*0.5+0.5

# populate panningtables
b1 = np.zeros((1024,))
b2 = np.zeros((1024,))
b3 = np.zeros((1024,))

b1[:512] = sine_gain
b2[:512] = sine_gain
b3[:512] = sine_gain

b2[512:] = b2[:512]*0.5
b3[512:] = b3[:512]

# input frequency
f = np.sin(np.linspace(0, 2*np.pi*f, fs))

# rotations (circular panning) using each buffer type
r1 = np.interp(np.linspace(0,fr,fs)%1., np.linspace(0,1,1024), b1)
r2 = np.interp(np.linspace(0,fr,fs)%1., np.linspace(0,1,1024), b2)
r3 = np.interp(np.linspace(0,fr,fs)%1., np.linspace(0,1,1024), b3)


fig, ax = plt.subplots(4,1)
[ax[0].plot(np.linspace(0, 2*np.pi, len(i[0])), i[0], linestyle=i[1]) for i in zip((b1, b2, b3), ('dashed', '-.', ':'))]
ax[0].set_xticks(np.linspace(0, np.pi*2, 7), ["0", r'$\pi/3$', r'$2\pi/3$', r'$\pi$', r'$4\pi/3$', r'$5\pi/3$', r'$2\pi$'])
ax[0].set_yticks([0, 0.5, 1], ["0.0", "0.5", "1.0"])
ax[0].set_xlabel('Virtual source position')
ax[0].set_ylabel('Gain')
ax[0].grid(True, which='both')
ax[0].set_xticks(np.linspace(0, np.pi*2, 7)[:-1]+np.pi*1/6, minor=True)
[ax[0].axvline(np.pi*i/2, linestyle=":", color='grey', lw=2) for i in range(4)]
ax[0].axvline(np.pi*1/2, linestyle="-", color='black', lw=2)
[i[0].magnitude_spectrum(i[1], Fs=fs) for i in zip(ax[1:], (f*r1, f*r2, f*r3))]
[i.set_xticks(np.linspace(0, 24000, 9),[]) for i in ax[1:]]
[i.set_xlabel('') for i in ax[1:]]
ax[3].set_xticks(np.linspace(0, 24000, 9), np.linspace(0, 24, 9).astype('int'))
ax[3].set_xlabel('Frequency (Hz)')
[i.set_ylim([0,0.26]) for i in ax[1:]]
[i.set_yticks([0, 0.1, 0.2]) for i in ax[1:]]
[i.set_ylabel('Magnitude') for i in ax[1:]]
linetypes = [('C0', 'dashed'), ('C1', '-.'), ('C2', ':')]
for i, e in enumerate(linetypes):
	dashed = Line2D([0], [0], color=e[0], linestyle=e[1])
	handles, labels = ax[i+1].get_legend_handles_labels()
	handles.extend([dashed])
	ax[i+1].legend(handles=handles)
[i.grid(True, which='both') for i in ax[1:]]
[i.set_xticks(np.linspace(1000, 23000, 23)[np.mod(np.arange(1,24), 3)!=0], minor=True) for i in ax[1:]]
plt.tight_layout()
plt.show()




