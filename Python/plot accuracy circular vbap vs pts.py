import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, optimize
from datetime import datetime
from scipy.interpolate import make_interp_spline

# Panning table synthesis (uniform distribution)
def pts(fs=48000, length=1, bits=512, channels=3, hori_f=1, interpolation='none'):

    def np_vbap(channels=3, azi=0.): # numpy version of VBAP
        azi = azi if isinstance(azi, np.ndarray) else np.array([azi]) 

        # create VBAP speaker circle
        speaker_array = np.linspace(0, np.pi*2, channels+1)[:-1]

        # compute distance to each speaker in each azimuth position
        azi_array = np.array([azi,] * channels) % (np.pi*2)
        distances = np.abs(speaker_array.reshape((channels, 1)) - azi_array)
        distance_circular = np.pi-np.abs(np.pi-distances)

        # find active speakers
        speaker_index = np.arange(channels)
        active_indices = [[asl for _, asl in sorted(zip(dc, speaker_index))][:2] for dc in distance_circular.T]

        # Helper function to normalize gain factors to garantie power conservation.
        def _normalize_gains(gains, vol_normalization=1): # vol_normalization between 0-1 
            return gains * np.sqrt(vol_normalization / np.sum(gains ** 2))

        azi_vector = np.asarray([-np.sin(azi), np.cos(azi)])

        gains = np.zeros((len(azi), channels))
        base_vectors = np.asarray([-np.sin(speaker_array), np.cos(speaker_array)])
        for i, active_index in enumerate(active_indices):
            active_base_vectors = base_vectors[:, active_index]
            inverted_base = np.linalg.inv(active_base_vectors)
            gains[i, active_index] = _normalize_gains(inverted_base @ azi_vector[:, i])

        return gains

    # pts buffer
    x = np.linspace(0, 1, bits)
    y = np_vbap(channels=channels, azi=x*2*np.pi)
    
    if interpolation=='cubic':
    #     x = np.insert(x, 1, x[1]/2)
    #     x = np.insert(x, -1, x[-2]+x[1]) # we already inserted x[1]/2 at position 1 in the line above! So, don't divide by 2 again!
    #     y = np.insert(y, 1, np_vbap(channels=channels, azi=x[1]), axis=0)
    #     y = np.insert(y, -1, np_vbap(channels=channels, azi=x[-2]), axis=0)
    #     wt = make_interp_spline(x, y, bc_type='not-a-knot')
        wt = make_interp_spline(x, y, bc_type='not-a-knot')

    # use constant DC offset signal
    input_audio = np.ones(int(fs*length))
    output_audio = np.zeros((channels, int(fs*length)))

    start=datetime.now()

    # constant circular rotation
    wr_h = 2 * np.pi * np.linspace(0, length, int(fs*length)) * hori_f

    # compute PTS
    sawtooth = signal.sawtooth(wr_h) * 0.5 + 0.5
    for j in range(channels):
        if interpolation=='none':
            horizontal_panning = y[(sawtooth*(bits-1)).astype(int), j]
        elif interpolation=='linear':
            horizontal_panning = np.interp(sawtooth, x, y[:,j])
        elif interpolation=='cubic':
            horizontal_panning = wt(sawtooth)[:, j]

        output_audio[j][:] = input_audio * horizontal_panning

    return output_audio, datetime.now()-start

# audio rate VBAP (uniform distribution)
def vbap_ar(fs=48000, length=1., channels=3, hori_f=1):
    x = np.linspace(0, length, int(fs*length))
    ch = np.sum(np.array(channels))

    # create VBAP speaker circle
    speaker_array = np.linspace(0, np.pi*2, channels+1)[:-1]

    # Helper function to normalize gain factors to garantie power conservation.
    def _normalize_gains(gains, vol_normalization=1): # vol_normalization between 0-1 
        return gains * np.sqrt(vol_normalization / np.sum(gains ** 2))

    # use constant DC offset signal
    input_audio = np.ones(int(fs*length))
    output_gains = np.zeros((int(fs*length), channels))

    start=datetime.now()

    # constant circular rotation
    wr_h = 2 * np.pi * np.linspace(0, length, int(fs*length)) * hori_f

    # to compute VBAP, we need the panning function as a continuous vector (array of vectors)
    wr_h_vector = np.asarray([-np.sin(wr_h), np.cos(wr_h)])

    # compute distance to each speaker in each audio sample
    wr_h_array = np.array([wr_h,]*channels) % (np.pi*2)
    distances = np.abs(speaker_array.reshape((channels, 1)) - wr_h_array)
    distance_circular = np.pi-np.abs(np.pi-distances)

    # find active speakers (2 in each sample in 2D)
    speaker_index = np.arange(channels)
    active_indices = [[asl for _, asl in sorted(zip(dc, speaker_index))][:2] for dc in distance_circular.T]

    # compute VBAP
    base_vectors = np.asarray([-np.sin(speaker_array), np.cos(speaker_array)])
    for i, active_index in enumerate(active_indices):
        active_base_vectors = base_vectors[:, active_index]
        inverted_base = np.linalg.inv(active_base_vectors)
        output_gains[i, active_index] = _normalize_gains(inverted_base @ wr_h_vector[:, i])

    output_audio = (output_gains * np.repeat(input_audio[:, np.newaxis], channels, axis=1)).T

    return output_audio, datetime.now()-start


# example usage
if __name__ == "__main__":
    channels = 3

    vbap_times = []
    print("--- Computing VBAP: (computation times in Python should not be considered as a valid comparison!)")
    for i in range(21):
        vbap_output, vbap_time = vbap_ar(channels=channels)
        print("VBAP done in: %f"%vbap_time.total_seconds())
        vbap_times.append(vbap_time.total_seconds())
    print("--- VABP average: %f"%np.mean(np.array(vbap_times)))

    print("\n--- Computing PTS: (computation times in Python should not be considered as a valid comparison!)")
    pts_times = []
    maxs = {}
    stds = {}
    avgs = {}
    x = np.arange(0,21)
    for interpolation in ['none', 'linear']:#, 'cubic']: # not-a-knot cubic interpolation not correctly implemented
        maxs[interpolation] = []
        stds[interpolation] = []
        avgs[interpolation] = []
        for i in x:
            print("Computing buffer size 2^i with i=%i: %i"%(i,2**i))
            pts_output, pts_time = pts(bits=2**i, channels=channels, interpolation=interpolation)
            print("  PTS done in: %f" %pts_time.total_seconds())
            pts_times.append(pts_time.total_seconds())

            difference = np.abs(pts_output - vbap_output)
            maxs[interpolation].append(np.max(difference))
            stds[interpolation].append(np.std(difference))
            avgs[interpolation].append(np.average(difference))

            # fig, ax = plt.subplots(2,1)
            # ax[0].plot(pts_output[0, :])
            # ax[0].plot(vbap_output[0, :])
            # ax[1].plot(difference[0,:])
            # fig.suptitle("Max error: %f"%maxs[interpolation][-1])
            # plt.show()
    print("--- PTS average: %f"%np.mean(np.array(pts_times)))

    # -- plot VBAP curves
    p = np.linspace(0, 2*np.pi, len(pts_output[0,:]))
    fig, ax = plt.subplots(channels,1)
    for i in range(channels):
        ax[i].plot(p, pts_output[i,:])
        ax[i].set_ylabel('Gain')
        ax[i].xaxis.set_tick_params(labelbottom=False)
        ax[i].set_xticks([0, np.pi*1/3, np.pi*2/3, np.pi, np.pi*4/3, np.pi*5/3, np.pi*2], ['$0$', '$\pi/3$', '$2\pi/3$', '$\pi$', '$4\pi/3$', '$5\pi/3$', '$2\pi$'])
    ax[-1].xaxis.set_tick_params(labelbottom=True)
    ax[-1].set_xlabel('Virtual source position')
    plt.show()

    # -- plot accuracy per power
    fig, ax1 = plt.subplots()
    ax1.plot(x, avgs['none'], '.', label='No interpolation')
    ax1.plot(x, avgs['linear'], 'x', label='Linear interpolation')
    # ax1.plot(x, avgs['cubic'], '+', label='Cubic')
    ax1.set_yscale('log')
    ax1.set_xlabel('Panning table size')
    ax1.set_ylabel('Avgerage error')
    ax1.set_xticks([0, 4, 8, 12, 16, 20], ['$2^0$', '$2^4$', '$2^8$', '$2^{12}$', '$2^{16}$', '$2^{20}$'])
    ax1.minorticks_on()
    ax1.grid(color='grey', linewidth=0.5, which='major')
    ax1.grid(color='grey', linewidth=0.3, which='minor')
    ax1.legend(loc="upper right")
    plt.show()


    # -- plot improvement of no vs. linear interpolation and show fitted curve
    fig, ax1 = plt.subplots()
    y = np.array(avgs['linear'])/np.array(avgs['none'])
    ax1.plot(x, y)
    # opt_res = optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  x,  y,  p0=(4, 0.1))[0] # logarithmic fit
    # y_fit = opt_res[0] * np.exp(opt_res[1] * x) # logarithmic fit didn't quite match the curve
    y_fit = 1.12165053 * np.exp(-0.66883177 * x) # slight modification for better fit
    ax1.plot(x, y_fit, '--')
    ax1.set_yscale('log')
    ax1.set_xlabel('Panning table size')
    ax1.set_ylabel('Avgerage error')
    ax1.set_ylabel('Ratio of error using no vs. linear interpolation')
    ax1.set_xticks([0, 4, 8, 12, 16, 20], ['$2^0$', '$2^4$', '$2^8$', '$2^{12}$', '$2^{16}$', '$2^{20}$'])
    ax1.minorticks_on()
    ax1.grid(color='grey', linewidth=0.5, which='major')
    ax1.grid(color='grey', linewidth=0.3, which='minor')
    ax1.set_title('Fitted curve at ~$1.12e^{-0.67x}$')
    plt.show()


