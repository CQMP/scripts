"""Utility functions for writing ALPScore-compatible HDF5 archives.

Simple library to create (the rather involved) metadata information as required
by ALPScore and outlined in the `ALPScore H5GF archive specification
<https://github.com/ALPSCore/H5GF/blob/master/H5GF.rst>` from outside of
ALPS.

To create the metadata for a quantity manually, you can follow this recipe::

    import h5py
    import numpy as np
    from h5archive import gf

    f = h5py.File('test.hdf5')
    qtty = gf.new_quantity(f, 'G0')
    gf.new_mesh_group(qtty, 3)
    gf.new_matsubara_mesh(qtty, 1, 30, beta=10.)
    gf.new_kspace_mesh(qtty, 2, [[0.]])
    gf.new_index_mesh(qtty, 3, 2)
    gf.new_inftail(f_g0iw, [np.ones((1,2))])
    gf.new_data(f_g0iw, np.zeros((30,1,2), complex))
"""
import numpy as _np
import h5py as _h5py

# TODO: make a nice interface for writing a quantity in one go.  Maybe with
#       objects and classes and all those things such as that.

def new_quantity(parent, name, creator=None):
    """Creates a new group for a stored quantity and returns a handle to it"""
    version = 0, 1
    ref_url = 'https://github.com/ALPSCore/H5GF/blob/master/H5GF.rst'
    if creator is None:
        creator = 'DontBeLazyAndChangeTheDefaults(TM) 0.0.0.1'

    qtty = parent.create_group(name)
    try:
        group = qtty.create_group("version")
        group.create_dataset("major", data=version[0], dtype="<i4")
        group.create_dataset("minor", data=version[1], dtype="<i4")
        group.create_dataset('originator', data=_make_string(creator))
        group.create_dataset('reference', data=_make_string(ref_url))
    except:
        # Remove newly created group in case something goes wrong in order to
        # provide 'strong' exception guarantee (being able to try again). Note
        # that this may break in cases of parallel writes (Parallel HDF5), but
        # is unlikely to because (1.) not more than one core typically writes
        # metadata and (2.) H5Py does not support parallel writes yet :)
        del parent['name']
        raise
    else:
        return qtty

def new_inftail(group, moments=None):
    """Creates a simple Matsubara frequency tail by storing moments

    As `moments`, expects a list of moments `A_n` for `n = 1, 2, ...`.
    """
    tail_type = 'INFINITY_TAIL'
    if moments is None:
        min_order = -1
        max_order = -1
        moments = []
    else:
        min_order = 1
        max_order = len(moments)

    tail = group.create_group("tail")
    try:
        tail.create_dataset('descriptor', data=_make_string(tail_type))
        tail.create_dataset('min_tail_order', data=min_order, dtype="<i4")
        tail.create_dataset('max_tail_order', data=max_order, dtype="<i4")
        for order, moment in enumerate(moments, start=min_order):
            tail.create_dataset(str(order), data=moment)
    except:
        del group['tail']
        raise
    else:
        return tail

def new_data(group, data):
    """Creates the data, handling complex data accurately"""
    data = _np.asarray(data)
    is_complex = _np.issubdtype(data.dtype, _np.complexfloating)

    if is_complex:
        data_write = _np.empty(data.shape + (2,), data.real.dtype)
        data_write[..., 0] = data.real
        data_write[..., 1] = data.imag
    else:
        data_write = data

    dset = group.create_dataset("data", data=data_write)
    try:
        if is_complex:
            dset.attrs["__complex__"] = 1
    except:
        del group['data']
        raise
    else:
        return dset

def new_mesh_group(group, naxes):
    """Creates the base group for the mesh information.

    This function must be called before calling any of the other mesh functions
    (`new_tau_mesh()`, `new_index_mesh()`, etc.)
    """
    mesh = group.create_group("mesh")
    try:
        mesh.create_dataset('N', data=naxes, dtype="<i4")
    except:
        del group['mesh']
        raise
    else:
        return mesh

def new_tau_mesh(group, axis, ntau, beta, mtype='mesh', stat='fermi'):
    """Creates a new mesh in imaginary time from [0, beta]"""
    kind = 'IMAGINARY_TIME'
    statistics = {'bose': 0, 'fermi': 1}[stat]
    if mtype == 'mesh':
        centres = False
        endpoint = True
        values = float(beta)/(ntau - 1) * _np.arange(0, ntau)
    elif mtype == 'bins':
        centres = True
        endpoint = False
        values = float(beta)/(2 * ntau) * _np.arange(1, 2*ntau, 2)
    else:
        raise ValueError("Unknown mesh type")

    mesh = group.create_group('mesh/%d' % axis)
    try:
        mesh.create_dataset('kind', data=_make_string(kind))
        mesh.create_dataset('N', data=ntau, dtype="<i4")
        mesh.create_dataset('statistics', data=statistics, dtype="<i4")
        mesh.create_dataset('beta', data=beta, dtype="<f8")
        mesh.create_dataset('last_point_included', data=endpoint, dtype="<i4")
        mesh.create_dataset('half_point_mesh', data=centres, dtype="<i4")
        mesh.create_dataset('points', data=values, dtype="<f8")
    except:
        del group['mesh/%d' % axis]
        raise
    else:
        return mesh

def new_matsubara_mesh(group, axis, niw, beta, full=False, stat='fermi'):
    """Creates a new mesh for Matsubara frequencies"""
    kind = 'MATSUBARA'
    full = bool(full)

    if stat == 'fermi':
        statistics = 1
        values = _np.pi/beta * _np.arange(-2*niw * full + 1, 2*niw, 2)
    elif stat == 'bose':
        statistics = 0
        values = _np.pi/beta * _np.arange(-2*(niw - 1) * full, 2*niw, 2)
    else:
        raise ValueError("Statistics must be either 'fermi' or 'bose'")

    mesh = group.create_group('mesh/%d' % axis)
    try:
        mesh.create_dataset('kind', data=_make_string(kind))
        mesh.create_dataset('N', data=niw, dtype="<i4")
        mesh.create_dataset('statistics', data=statistics, dtype="<i4")
        mesh.create_dataset('beta', data=beta, dtype="<f8")
        mesh.create_dataset('positive_only', data=not full, dtype="<i4")
        mesh.create_dataset('points', data=values, dtype="<f8")
    except:
        del group['mesh/%d' % axis]
        raise
    else:
        return mesh

def new_index_mesh(group, axis, nindex):
    """Creates a new mesh for a generic index such as spin"""
    kind = 'INDEX'

    mesh = group.create_group('mesh/%d' % axis)
    try:
        mesh.create_dataset('kind', data=_make_string(kind))
        mesh.create_dataset('N', data=nindex, dtype="<i4")
    except:
        del group['mesh/%d' % axis]
        raise
    else:
        return mesh

def new_kspace_mesh(group, axis, kpoints):
    """Creates a new mesh for a set of points in momentum space"""
    kind = 'MOMENTUM_INDEX'

    mesh = group.create_group('mesh/%d' % axis)
    try:
        mesh.create_dataset('kind', data=_make_string(kind))
        mesh.create_dataset('points', data=kpoints, dtype="<f8")
    except:
        del group['mesh/%d' % axis]
        raise
    else:
        return mesh

_hdf5_varstring_dtype = _h5py.special_dtype(vlen=str)

def _make_string(s):
    """Return object suitable for writing to HDF5 file from string"""
    return _np.array(s, dtype=_hdf5_varstring_dtype)
