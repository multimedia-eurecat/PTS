import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def vbap(phi, phi_0):
    return np.sqrt((np.tan(phi_0)**2 + 2*np.tan(phi_0)*np.tan(phi) + np.tan(phi)**2) / (2*(np.tan(phi_0)**2+np.tan(phi)**2)))

channels = 10
panningtable_size = 1024
phi_0 = np.pi/np.max(channels)

panningtables = np.zeros((channels, panningtable_size))

phi_0 = np.pi/np.max(channels)
panningtable_range = int(panningtable_size / channels)
for i in range(channels):
    position = ((panningtable_size-panningtable_range)+panningtable_range*i) % panningtable_size
    panningtables[i][position:position+panningtable_range] = vbap(np.linspace(-phi_0, phi_0, panningtable_range), phi_0)
    position = (position+panningtable_range) % panningtable_size
    panningtables[i][position:position+panningtable_range] = vbap(np.linspace(phi_0, -phi_0, panningtable_range), phi_0)

fig, ax = plt.subplots(channels,1) #, figsize=(1920/300,1043/300), dpi=300)
for i in range(channels):
    ax[i].plot(panningtables[i])
    ax[i].set_xticks([0,256,512,768,1024],[])
    ax[i].set_yticks([0, 0.5, 1], ['0.0', '0.5', '1.0'])
    ax[i].set_xlabel('')
    ax[i].set_ylabel('Gain')
    ax[i].grid(True, which='both')
    ax[i].set_xlim(0,1024)
    for j in range(channels):
        pos = j*panningtable_size/channels
        ax[i].axvline(1 if pos == 0 else pos, linestyle='-' if i==j else ':', color='black' if i==j else 'grey', lw=2)
ax[channels-1].set_xticks([0,256,512,768,1024],[0,128,256,384,512])
ax[channels-1].set_xlabel('Azimuth points')
# plt.tight_layout(rect=(0, 0, 1, 0.88))
plt.show()