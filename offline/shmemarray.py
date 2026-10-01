import mmap

import numpy
from packaging.version import Version

NUMPY_LT_2_0 = not (Version(numpy.__version__) < Version("2.0.dev"))


def empty_like(array, dtype=None):
    """Create a shared memory array from the shape of array."""
    array = numpy.asarray(array)
    if dtype is None:
        dtype = array.dtype
    return anonymousmemmap(array.shape, dtype)


def empty(shape, dtype="f8"):
    """Create an empty shared memory array."""
    return anonymousmemmap(shape, dtype)


def full_like(array, value, dtype=None):
    """Create a shared memory array with the same shape and type
    as a given array, filled with `value`."""
    shared = empty_like(array, dtype)
    shared[:] = value
    return shared


def full(shape, value, dtype="f8"):
    """Create a shared memory array of given shape and type,
    filled with `value`."""
    shared = empty(shape, dtype)
    shared[:] = value
    return shared


def copy(a):
    """Copy an array to the shared memory.

    Notes
    -----
    copy is not always necessary because the private memory is
    always copy-on-write.

    Use :code:`a = copy(a)` to immediately dereference the old 'a'
    on private memory
    """
    shared = anonymousmemmap(a.shape, dtype=a.dtype)
    shared[:] = a[:]
    return shared


def fromiter(iter, dtype, count=None):
    return copy(numpy.fromiter(iter, dtype, count))


try:
    # numpy >= 1.16
    _unpickle_ctypes_type = numpy.ctypeslib.as_ctypes_type(numpy.dtype("|u1"))
except Exception:
    # older version numpy < 1.16
    _unpickle_ctypes_type = numpy.ctypeslib._typecodes["|u1"]


def __unpickle__(ai, dtype):
    dtype = numpy.dtype(dtype)
    tp = _unpickle_ctypes_type

    # The following lb ub logic replaces ctypeslib.as_ctypes()
    # that does not support strides.
    if ai["strides"]:
        lb = 0
        ub = dtype.itemsize
        for s, t in zip(ai["strides"], ai["shape"]):
            if t == 0:  # there is no data if any shape is 0.
                lb, ub = 0, 0
                break
            if s < 0:
                lb = lb + (t - 1) * s
            else:
                ub = ub + (t - 1) * s
    else:
        lb = 0
        ub = dtype.itemsize
        for s in ai["shape"]:
            ub = ub * s

    # grab a flat char array at the sharemem address,
    # covering the memory region ai refers to.
    ra = (tp * (ub - lb)).from_address(ai["data"][0] + lb)

    # view it as what it should look like. here we assume numpy
    # will not do crazy things like modifying gaps due to striding.
    shm = numpy.ndarray(
        buffer=ra, dtype=dtype, offset=-lb,
        strides=ai["strides"], shape=ai["shape"]
    ).view(type=anonymousmemmap)
    return shm


class anonymousmemmap(numpy.memmap):
    """Arrays allocated on shared memory.

    The array is stored in an anonymous memory map
    that is shared between child-processes.

    """

    def __new__(subtype, shape, dtype=numpy.uint8, order="C"):

        descr = numpy.dtype(dtype)
        _dbytes = descr.itemsize

        shape = numpy.atleast_1d(shape)
        size = 1
        for k in shape:
            size *= k

        bytes = int(size * _dbytes)

        if bytes > 0:
            mm = mmap.mmap(-1, bytes)
        else:
            mm = numpy.empty(0, dtype=descr)
        self = numpy.ndarray.__new__(
            subtype, shape, dtype=descr, buffer=mm, order=order
        )
        self._mmap = mm
        return self

    def __array_wrap__(self, outarr, context=None, return_scalar=False):
        # after ufunc this won't be on shm!
        if NUMPY_LT_2_0:
            return numpy.ndarray.__array_wrap__(
                self.view(numpy.ndarray), outarr, context
            )
        else:
            return numpy.ndarray.__array_wrap__(
                self.view(numpy.ndarray), outarr, context, return_scalar
            )

    def __reduce__(self):
        return __unpickle__, (self.__array_interface__, self.dtype)

