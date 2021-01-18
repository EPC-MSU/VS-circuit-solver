from ctypes import CDLL, Structure, Array, c_ubyte, \
    c_double, c_size_t, POINTER
from platform import system
import numpy as np
VOLTAGE_AMPL = 12.
R_CS = 475.
CURRENT_AMPL = (VOLTAGE_AMPL / R_CS * 1000)


def _get_dll():
    if system() == "Linux":
        return CDLL("libivcmp.so")
    elif system() == "Windows":
        return CDLL("ivcmp.dll")
    else:
        raise NotImplementedError("Unsupported platform {0}".format(system()))


lib = _get_dll()
MAX_NUM_POINTS = 1000


class _IterableStructure(Structure):
    def __iter__(self):
        return (getattr(self, n) for n, t in self._fields_)


def _normalize_arg(value, desired_ctype):
    from collections import Sequence

    if isinstance(value, desired_ctype):
        return value
    elif issubclass(desired_ctype, Array) and isinstance(value, Sequence):
        member_type = desired_ctype._type_

        if desired_ctype._length_ < len(value):
            raise ValueError()

        if issubclass(member_type, c_ubyte) and isinstance(value, bytes):
            return desired_ctype.from_buffer_copy(value)
        elif issubclass(member_type, c_ubyte) and isinstance(value, bytearray):
            return value
        else:
            return desired_ctype(*value)
    else:
        return value


class IvCurve(_IterableStructure):
    _fields_ = (
        ("voltages", c_double*MAX_NUM_POINTS),
        ("currents", c_double*MAX_NUM_POINTS),
        ("length", c_size_t)
    )

    def __init__(self):
        self.length = MAX_NUM_POINTS


def SetMinVC(min_var_v, min_var_c):
    lib_func = lib.SetMinVC
    lib_func.argtype = c_double, c_double
    lib_func(c_double(min_var_v),  c_double(min_var_c))


def CompareIvc(first_iv_curve, second_iv_curve):
    lib_func = lib.CompareIVC
    lib_func.argtype = POINTER(c_double), POINTER(c_double), c_size_t, POINTER(c_double), POINTER(c_double), c_size_t
    lib_func.restype = c_double
    res = lib_func(first_iv_curve.voltages, first_iv_curve.currents, first_iv_curve.length,
                   second_iv_curve.voltages, second_iv_curve.currents, second_iv_curve.length)
    res = _normalize_arg(res, c_double)
    return res


if __name__ == "__main__":
    iv_curve = IvCurve()
    ivc_curve = IvCurve()
    for i in range(MAX_NUM_POINTS):
        iv_curve.voltages[i] = 0.47 * VOLTAGE_AMPL * np.sin(2 * 3.14 * i / MAX_NUM_POINTS)
        iv_curve.currents[i] = 0.63 * CURRENT_AMPL * np.sin(2 * 3.14 * i / MAX_NUM_POINTS)
        ivc_curve.voltages[i] = VOLTAGE_AMPL * np.sin(2 * 3.14 * i / MAX_NUM_POINTS)
        ivc_curve.currents[i] = CURRENT_AMPL * np.cos(2 * 3.14 * i / MAX_NUM_POINTS)
    SetMinVC(0, 0)
    f = CompareIvc(iv_curve, ivc_curve)
    print(f)
    # for i in range(MAX_NUM_POINTS):
    #     print(iv_curve.currents[i], iv_curve.voltages[i])
