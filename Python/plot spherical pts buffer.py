import numpy as np
import scipy.spatial as scyspat
import matplotlib.pyplot as plt

# VBAP code extracted from spaudiopy repository: https://github.com/chris-hld/spaudiopy

def deg2rad(deg):
    """Convert from degree [0, 360) to radiant [0, 2*pi)."""
    return deg % 360 / 180 * np.pi


def rad2deg(rad):
    """Convert from radiant [0, 2*pi) to degree [0, 360)."""
    return rad / np.pi * 180 % 360

def cart2sph(x, y, z, steady_colat=False):
    """Vectorized conversion of cartesian to spherical coordinates."""
    x = asarray_1d(x)
    y = asarray_1d(y)
    z = asarray_1d(z)
    r = np.sqrt(x**2 + y**2 + z**2)
    azi = np.arctan2(y, x)
    colat = np.arccos(z / r) if not steady_colat else \
            np.arccos(z / np.clip(r, 10e-15, None))
    return azi, colat, r

def sph2cart(azi, colat, r=1):
    """Vectorized conversion of spherical to cartesian coordinates."""
    azi = asarray_1d(azi)
    colat = asarray_1d(colat)
    r = asarray_1d(r)
    x = r * np.cos(azi) * np.sin(colat)
    y = r * np.sin(azi) * np.sin(colat)
    z = r * np.cos(colat)
    return x, y, z

def area_triangle(p1, p2, p3):
    """calculate area of any triangle given coordinates of its corners p."""
    return 0.5 * np.linalg.norm(np.cross((p2 - p1), (p3 - p1)))

def asarray_1d(a, **kwargs):
    """Squeeze the input and check if the result is one-dimensional.
    Returns *a* converted to a `numpy.ndarray` and stripped of
    all singleton dimensions.  Scalars are "upgraded" to 1D arrays.
    The result must have exactly one dimension.
    If not, an error is raised.
    """
    result = np.squeeze(np.asarray(a, **kwargs))
    if result.ndim == 0:
        result = result.reshape((1,))
    elif result.ndim > 1:
        raise ValueError("array must be one-dimensional")
    return result

class LoudspeakerSetup:
    """Creates a 'hull' object containing all information for further decoding.
    .. plot::
        :context: close-figs
        ls_setup = spa.decoder.LoudspeakerSetup(ls_x, ls_y, ls_z)
        ls_setup.pop_triangles(normal_limit=85, aperture_limit=90,
                               opening_limit=150)
        ls_setup.show()
    """

    def __init__(self, x, y, z, listener_position=None):
        """
        Parameters
        ----------
        x : array_like
        y : array_like
        z : array_like
        listener_position : (3,), cartesian, optional
            Offset, will be substracted from the loudspeaker positions.
        """
        self.x = asarray_1d(x)
        self.y = asarray_1d(y)
        self.z = asarray_1d(z)
        if listener_position is None:
            listener_position = [0, 0, 0]
        self.listener_position = asarray_1d(listener_position)

        # Listener position as origin
        self.x -= listener_position[0]
        self.y -= listener_position[1]
        self.z -= listener_position[2]
        # TODO: Better handling of this, e.g. not effective when updating hull:
        self.listener_position -= self.listener_position
        _, _, self.d = cart2sph(self.x, self.y, self.z)

        # amplitude decay exponent
        self.a = 1

        # Loudspeaker static gains
        self.ls_gains = np.ones(len(self.d))  # Available, but not used, yet

        # Triangulation of points
        self._hull = get_hull(self.x, self.y, self.z)
        self.points = self._hull.points
        self.npoints = self._hull.npoints
        self.nsimplex = self._hull.nsimplex
        self.vertices = self._hull.vertices
        self.simplices = self._hull.simplices
        self.simplices = sort_vertices(self.simplices)
        self.centroids = calculate_centroids(self)
        self.face_areas = calculate_face_areas(self)
        self.face_normals = calculate_face_normals(self)
        self.vertex_normals = calculate_vertex_normals(self)
        self.barycenter = calculate_barycenter(self)

        # All simplices enclosing listener valid per default for rendering
        self._encloses_listener = check_listener_inside(self,
                                                        self.listener_position)
        self.valid_simplices = self._encloses_listener

        # VBAP
        self.inverted_vbase = None  # populated by VBAP

        # see 'ambisonics_setup()'
        self.ambisonics_hull = []
        self.kernel_hull = []
        self.characteristic_order = None

        # some checks
        assert(len(self.d) == self.npoints)

    @classmethod
    def from_sph(cls, azi, colat, r=1, listener_position=None):
        """ Alternative constructor, using spherical coordinates in rad.
        Parameters
        ----------
        azi : array_like, spherical
        colat : array_like, spherical
        r : array_like, spherical
        listener_position : (azi, colat, r), spherical, optional
            Offset, will be substracted from the loudspeaker positions.
        """
        x, y, z = sph2cart(azi, colat, r)
        if listener_position is None:
            listener_position = [0, 0, 0]
        listener_position = asarray_1d(listener_position)
        listener_position = sph2cart(*listener_position)
        return cls(x, y, z, listener_position=listener_position)

    def is_simplex_valid(self, simplex):
        """Tests if simplex is in valid simplices (independent of orientation).
        """
        # find face in all faces
        in_s = np.isin(self.valid_simplices, simplex).sum(axis=-1) == 3
        return np.any(in_s)

    def pop_triangles(self, normal_limit=85, aperture_limit=None,
                      opening_limit=None, blacklist=None):
        """Refine triangulation by removing them from valid simplices.
        Bypass by passing 'None'.
        Parameters
        ----------
        normal_limit : float, optional
        aperture_limit : float, optional
        opening_limit : float, optional
        blacklist : list, optional
        """
        if normal_limit is not None:
            self.valid_simplices = check_normals(self, normal_limit)
        if aperture_limit is not None:
            self.valid_simplices = check_aperture(self, aperture_limit)
        if opening_limit is not None:
            self.valid_simplices = check_opening(self, opening_limit)
        if blacklist is not None:
            self.valid_simplices = apply_blacklist(self, blacklist)

    def get_characteristic_order(self):
        """Characteristic Ambisonics order."""
        if self.characteristic_order is None:
            N_e = characteristic_ambisonic_order(self)
            if N_e < 1:
                raise ValueError
        else:
            N_e = self.characteristic_order
        return N_e

    def ambisonics_setup(self, N_kernel=50, update_hull=False,
                         imaginary_ls=None):
        """Prepare loudspeaker hull for ambisonic rendering.
        Sets the `kernel_hull` as an n-design of twice `N_kernel`,
        and updates the ambisonic hull with an additional imaginary
        loudspeaker, if desired.
        Parameters
        ----------
        N_kernel : int, optional
        update_hull : bool, optional
        imaginary_ls : (L, 3), cartesian, optional
            Imaginary loudspeaker positions, if set to 'None' calls
            'find_imaginary_loudspeaker()' for 'update_hull'.
        Examples
        --------
        .. plot::
            :context: close-figs
            ls_setup.ambisonics_setup(update_hull=True)
            N_e = ls_setup.characteristic_order
            ls_setup.ambisonics_hull.show(title=f"Ambisonic Hull, $N_e={N_e}$")
        """
        self.characteristic_order = self.get_characteristic_order()
        if N_kernel is None:
            warn('Setting setup kernel order =', self.characteristic_order)
            N_kernel = self.characteristic_order
        if(not update_hull and imaginary_ls is not None):
            warn('Not updating hull but imaginary_ls position given.')

        ambi_ls = self.points
        if update_hull:
            if imaginary_ls is None:
                new_imaginary_ls = find_imaginary_loudspeaker(self)
                # add imaginary speaker to hull
                ambi_ls = np.vstack([ambi_ls, new_imaginary_ls])
                # mark imaginary speaker (last one)
                imaginary_ls_idx = ambi_ls.shape[0] - 1
            else:
                imaginary_ls = np.atleast_2d(imaginary_ls)
                assert(imaginary_ls.shape[1] == 3)
                # add imaginary loudspeaker(s) to hull
                ambi_ls = np.vstack([ambi_ls, imaginary_ls])
                # mark imaginary speaker (last one(s))
                imaginary_ls_idx = np.arange(ambi_ls.shape[0] -
                                             imaginary_ls.shape[0],
                                             ambi_ls.shape[0])
        else:
            imaginary_ls_idx = None

        # This new triangulation is now the rendering setup
        ambisonics_hull = LoudspeakerSetup(ambi_ls[:, 0],
                                           ambi_ls[:, 1],
                                           ambi_ls[:, 2])
        # mark imaginary speaker index
        ambisonics_hull.imaginary_ls_idx = imaginary_ls_idx
        # discretization hull
        virtual_speakers = grids.load_n_design(2 * N_kernel)
        # Avoid any extra calculation on this dense grid
        kernel_hull = get_hull(virtual_speakers[:, 0],
                               virtual_speakers[:, 1],
                               virtual_speakers[:, 2])

        del ambisonics_hull.ambisonics_hull
        del ambisonics_hull.kernel_hull
        self.ambisonics_hull = ambisonics_hull
        self.kernel_hull = kernel_hull
        self.kernel_hull.N_kernel = N_kernel

    def binauralize(self, ls_signals, fs, orientation=(0, 0), hrirs=None):
        """Create binaural signals that the loudspeaker signals produce on this
        setup (no delays).
        Parameters
        ----------
        ls_signals : (L, S) np.ndarray
            Loudspeaker signals.
        fs : int
        orientation : (azi, colat) tuple, optional
            Listener orientation offset (azimuth, colatitude) in rad.
        hrirs : sig.HRIRs, optional
        Returns
        -------
        l_sig : array_like
        r_sig : array_like
        """
        if hrirs is None:
            hrirs = io.load_hrirs(fs)
        assert(hrirs.fs == fs)
        ls_signals = np.atleast_2d(ls_signals)
        assert ls_signals.shape[0] == self.npoints, \
            'Provide signal per loudspeaker!'
        # distance attenuation
        relative_position = self.points - \
                            self.listener_position
        ls_azi, ls_colat, ls_r = cart2sph(*relative_position.T)
        ls_signals = np.diag(1 / ls_r ** self.a) @ ls_signals
        # convolve with hrir
        l_sig = np.zeros(ls_signals.shape[1] + len(hrirs) - 1)
        r_sig = np.zeros_like(l_sig)
        for ch, ls_sig in enumerate(ls_signals):
            if any(abs(ls_sig) > 10e-6):  # Gate at -100dB
                hrir_l, hrir_r = hrirs.nearest_hrirs(ls_azi[ch] -
                                                     orientation[0],
                                                     ls_colat[ch] -
                                                     orientation[1])
                # sum all loudspeakers
                l_sig += signal.convolve(ls_sig, hrir_l)
                r_sig += signal.convolve(ls_sig, hrir_r)
        return l_sig, r_sig

    def loudspeaker_signals(self, ls_gains, sig_in=None):
        """Render loudspeaker signals.
        Parameters
        ----------
        ls_gains : (S, L) np.ndarray
        sig_in : (S,) array like, optional
        Returns
        -------
        sig_out : (L, S) np.ndarray
        """
        ls_gains = np.atleast_2d(ls_gains)
        if sig_in is None:
            sig_in = np.ones(ls_gains.shape[0])
        sig_in = asarray_1d(sig_in)
        assert(ls_gains.shape[1] == len(self.points)), \
            'Provide gain per speaker!'
        return (sig_in[:, np.newaxis] * ls_gains).T

    def show(self, title='Loudspeaker Setup', **kwargs):
        """Plot hull object, calls plot.hull()."""
        plot.hull(self, title=title, **kwargs)


def get_hull(x, y, z):
    """Wrapper for scipy.spatial.ConvexHull."""
    return scyspat.ConvexHull(np.c_[x, y, z], incremental=False)

def sort_vertices(simplices):
    """Start the simplices with smallest vertex entry."""
    out = np.zeros_like(simplices)
    for i, face in enumerate(simplices):
        face = face[::-1]
        out[i, :] = np.roll(face, -np.argmin(face))
    return out

def calculate_centroids(hull):
    """Calculate centroid for each simplex."""
    centroids = np.zeros((len(hull.simplices), 3))
    for face_i, face in enumerate(hull.simplices):
        # extract vertices face
        v = hull.points[face, :]
        centroids[face_i, :] = np.mean(v, axis=0)
    return centroids

def calculate_face_areas(hull):
    """Calculate area for each simplex."""
    face_areas = np.zeros(len(hull.simplices))
    for face_i, face in enumerate(hull.simplices):
        v = hull.points[face, :]
        face_areas[face_i] = area_triangle(v[0, :], v[1, :], v[2, :])
    return face_areas


def calculate_face_normals(hull, eps=10e-6, normalize=False):
    """Calculate outwards pointing normal for each simplex."""
    face_normals = np.zeros((len(hull.simplices), 3))
    barycenter = np.mean(hull.points, axis=0)
    for face_i, face in enumerate(hull.simplices):
        # extract vertices face
        v = hull.points[face, :]
        centroid = np.mean(v, axis=0)
        # normal vector is cross product of two sides, initial point v0
        v_n = np.cross(v[1, :] - v[0, :], v[2, :] - v[0, :])
        # compare of face normal points in direction of barycenter
        criterion = np.dot(centroid + v_n - centroid, barycenter - centroid)
        # Make normal vector pointing outwards
        if criterion > eps:
            v_n = -v_n
        if normalize:
            v_n = v_n / np.linalg.norm(v_n)
        face_normals[face_i, :] = v_n
    return face_normals


def calculate_vertex_normals(hull, normalize=False):
    """Calculate normal for each vertex from simplices normals."""
    faces = hull.simplices
    vertex_normals = np.zeros([hull.npoints, 3])
    for p in hull.vertices:
        is_in_face = [p in row for row in faces]
        a = hull.face_areas[is_in_face]
        N = hull.face_normals[is_in_face]
        # weighted sum
        N_w = a[:, np.newaxis] * N
        vertex_n = np.sum(N_w, axis=0)
        if normalize:
            vertex_n = vertex_n / np.linalg.norm(vertex_n)
        vertex_normals[p, :] = vertex_n
    return vertex_normals


def calculate_barycenter(hull):
    """Barycenter of hull object."""
    return np.mean(hull.points, axis=0)

def check_listener_inside(hull, listener_position=None):
    """Return valid simplices for which the listener is inside the hull."""
    if listener_position is None:
        listener_position = hull.listener_position
    listener_position = np.asarray(listener_position)
    valid_simplices = []
    for face, centroid in zip(hull.simplices, hull.centroids):
        # centroid to listener
        v1 = listener_position - centroid
        # centroid to barycenter
        v2 = hull.barycenter - centroid
        # listener inside if both point in the same direction
        if np.dot(v1, v2) < 0:
            print("Listener not inside:", face)
        else:
            valid_simplices.append(face)
    return np.array(valid_simplices)

def vbap(src, hull, norm=2, valid_simplices=None, retain_outside=False,
         jobs_count=1):
    """Loudspeaker gains for Vector Base Amplitude Panning decoding.
    Parameters
    ----------
    src : (n, 3) numpy.ndarray
        Cartesian coordinates of n sources to be rendered.
    hull : LoudspeakerSetup
    norm : non-zero int, float
        Gain normalization norm, e.g. 1: anechoic, 2: reverberant
    valid_simplices : (nsimplex, 3) numpy.ndarray
        Valid simplices employed for rendering, defaults hull.valid_simplices.
    retain_outside : bool, optional
        Render on the 'ambisonic hull' to fade out amplitude.
    jobs_count : int or None, optional
        Number of parallel jobs, 'None' employs 'cpu_count'.
    Returns
    -------
    gains : (n, L) numpy.ndarray
        Panning gains for L loudspeakers to render n sources.
    References
    ----------
    Pulkki, V. (1997). Virtual Sound Source Positioning Using Vector Base
    Amplitude Panning. AES, 144(5), 357–360.
    Examples
    --------
    .. plot::
        :context: close-figs
        ls_setup = spa.decoder.LoudspeakerSetup(ls_x, ls_y, ls_z)
        ls_setup.pop_triangles(normal_limit=85, aperture_limit=90,
                               opening_limit=150)
        spa.plot.decoder_performance(ls_setup, 'VBAP')
        ls_setup.ambisonics_setup(update_hull=True)
        spa.plot.decoder_performance(ls_setup, 'VBAP', retain_outside=True)
        plt.suptitle('VBAP with imaginary loudspeaker')
    """
    if jobs_count is None:
        jobs_count = multiprocessing.cpu_count()
    if retain_outside:
        assert(valid_simplices is None)
        if hull.ambisonics_hull:
            hull = hull.ambisonics_hull
            if hull.imaginary_ls_idx is None:
                raise ValueError('No imaginary loudspeaker. Update hull!')
        else:
            raise ValueError('Run LoudspeakerSetup.ambisonics_setup() first!')

    if valid_simplices is None:
        valid_simplices = hull.valid_simplices

    # inverted LS vector base
    if hull.inverted_vbase is None:
        inverted_vbase = _invert_triplets(valid_simplices, hull.points)
        hull.inverted_vbase = inverted_vbase
    else:
        inverted_vbase = hull.inverted_vbase

    src = np.atleast_2d(src)
    assert(src.shape[1] == 3)
    src_count = src.shape[0]

    ls_count = hull.npoints

    gains = np.zeros([src_count, ls_count])

    if (jobs_count == 1) or (src_count < 10):
        for src_idx in range(src_count):
            for face_idx, ls_base in enumerate(inverted_vbase):
                # projecting src onto loudspeakers
                projection = np.dot(ls_base, src[src_idx, :])
                # normalization
                projection /= np.linalg.norm(projection, ord=norm)
                if np.all(projection > -10e-6):
                    assert(np.count_nonzero(projection) <= 3)
                    # print(f"Source {src_idx}: Gains {projection}")
                    gains[src_idx, valid_simplices[face_idx]] = projection
                    break  # found valid gains
    else:
        logging.info("Using %i processes..." % jobs_count)
        # preparation
        shared_array_shape = np.shape(gains)
        _arr_base = _create_shared_array(shared_array_shape)
        _arg_itr = zip(range(src_count),
                       repeat(src), repeat(inverted_vbase),
                       repeat(valid_simplices), repeat(norm))
        # execute
        with multiprocessing.Pool(processes=jobs_count,
                                  initializer=_init_shared_array,
                                  initargs=(_arr_base,
                                            shared_array_shape,)) as pool:
            pool.starmap(_vbap_gains_single_source, _arg_itr)
        # reshape
        gains = np.frombuffer(_arr_base.get_obj()).reshape(
                                shared_array_shape)

    # Distance compensation
    gains = (hull.d[np.newaxis, :] ** hull.a) * gains
    if retain_outside:
        # remove imaginary loudspeaker
        gains = np.delete(gains, hull.imaginary_ls_idx, axis=1)
    return gains

def _invert_triplets(simplices, points):
    """Invert loudspeaker triplets."""
    inverted_ls_triplets = []
    for face in simplices:
        # extract vertices face (valid LS positions)
        v = points[face, :]
        v_inv = np.linalg.lstsq(v, np.eye(3), rcond=None)[0]
        inverted_ls_triplets.append(v_inv.T)
    return inverted_ls_triplets

def vogel_points(N):
    golden_angle = np.pi * (3 - np.sqrt(5))
    phi = golden_angle * np.arange(N)
    z = np.linspace(1-1.0/N, 1.0/N-1, N)
    radius = np.sqrt(1-z*z)
    
    pts = np.zeros((N, 3))
    pts[:,0] = radius * np.cos(phi)
    pts[:,1] = radius * np.sin(phi)
    pts[:,2] = z
    
    return pts

# example usage
if __name__ == "__main__":
    # -- Layout: simple 5 channel setup: Quad + VoG
    # ls_dirs = np.array([[-45, 45, 180-45, 180+45, 0],[0,0,0,0, 90]])
    # ls_dirs[1, :] = 90 - ls_dirs[1, :]
    # ls_x, ls_y, ls_z = sph2cart(deg2rad(ls_dirs[0, :]), deg2rad(ls_dirs[1, :]))

    # -- Layout: Vogel distribution
    #            Vogel allows wrapping a single scalar number of loudspeakers in an even distribution on the sphere
    # vogel = vogel_points(30)
    # ls_x, ls_y, ls_z = vogel[:,0], vogel[:,1], vogel[:,1]

    # -- Layout: Hamasaki 22.2 setup
    ls_dirs = np.array([[0, 0, 45, 90, 110, 180, -110, -90, -45, 0, 23, 45, 90, 110, 180, -110, -90, -45, -23, 0, 45, -45],
                        [90, 40, 40, 40, 40, 40, 40, 40, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -25, -25, -25]])
    ls_dirs[1, :] = 90 - ls_dirs[1, :]
    ls_x, ls_y, ls_z = sph2cart(deg2rad(ls_dirs[0, :]), deg2rad(ls_dirs[1, :]))

    # -- produce panningtable
    ls_setup = LoudspeakerSetup(ls_x, ls_y, ls_z, [0, 0, 0])
    bits = 512
    channel = 2 # which channel to visualise
    pt = np.zeros((bits//2, bits)) # 2D panningtable
    for el in range(bits//2):
        for az in range(bits):
            pos = np.array(sph2cart(az/bits*2*np.pi, el/(bits//2+1)*np.pi))
            print("  Measuring point on sphere x: %7.4f - y: %7.4f - z: %7.4f"%(pos[0], pos[1], pos[2]), end="\r")
            pt[el, az] = vbap(pos.T, ls_setup)[0, channel]

    # -- plot spherical gains
    ip = plt.imshow(pt, cmap='Greys')
    for ps in ls_dirs.T[1:]:
        plt.scatter((ps[0]/360+1)%1 * bits, ps[1]/360 * bits, color='black')
    plt.scatter((ls_dirs.T[channel][0]/360+1)%1 * bits, ls_dirs.T[channel][1]/360 * bits, facecolors='none', edgecolors='white')
    plt.xlim(0, bits)
    plt.ylim(bits//2, 0)
    plt.xticks(ticks=[0, bits//4, bits//2, bits*3//4, bits])
    plt.yticks(ticks=[0, bits//4, bits//2])
    plt.xlabel('Azimuth points')
    plt.ylabel('Elevation points')
    plt.minorticks_on()
    plt.grid(visible=True, which='both', axis='both')
    ax = plt.gca()
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(ip, cax=cax)
    cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    cbar.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    cbar.ax.set_ylabel('Gain', rotation=270)
    cbar.ax.get_yaxis().labelpad = 11
    plt.show()






