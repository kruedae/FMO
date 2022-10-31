#!python
#cython: language_level=3
# This file is generated automatically by QuTiP.
# (C) 2011 and later, QuSTaR
import numpy as np
cimport numpy as np
cimport cython
np.import_array()
cdef extern from "numpy/arrayobject.h" nogil:
    void PyDataMem_NEW_ZEROED(size_t size, size_t elsize)
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
    void PyDataMem_FREE(void * ptr)
from qutip.cy.interpolate cimport interp, zinterp
from qutip.cy.math cimport erf, zerf
cdef double pi = 3.14159265358979323
from qutip.cy.brtools cimport (dense_add_mult, ZHEEVR, dense_to_eigbasis,
        vec_to_eigbasis, vec_to_fockbasis, skew_and_dwmin,
        diag_liou_mult, spec_func, farray_alloc)
from qutip.cy.brtools cimport (cop_super_mult, br_term_mult)
include '/home/kennet/anaconda3/envs/bofin_env/lib/python3.8/site-packages/qutip/cy/complex_math.pxi'

cdef complex spectral0(double w, double t): return  2*pi* 2.0 * 6597344572538.565 / (pi * 6024096385542.169 * 2.5448096943108576e-14)  if (w==0) else 2*pi*(2.0*6597344572538.565*6024096385542.169 *w /(pi*(w**2+6024096385542.169**2))) * ((1/(exp((w) * 2.5448096943108576e-14)-1))+1)
cdef complex spectral1(double w, double t): return  2*pi* 2.0 * 6597344572538.565 / (pi * 6024096385542.169 * 2.5448096943108576e-14)  if (w==0) else 2*pi*(2.0*6597344572538.565*6024096385542.169 *w /(pi*(w**2+6024096385542.169**2))) * ((1/(exp((w) * 2.5448096943108576e-14)-1))+1)
cdef complex spectral2(double w, double t): return  2*pi* 2.0 * 6597344572538.565 / (pi * 6024096385542.169 * 2.5448096943108576e-14)  if (w==0) else 2*pi*(2.0*6597344572538.565*6024096385542.169 *w /(pi*(w**2+6024096385542.169**2))) * ((1/(exp((w) * 2.5448096943108576e-14)-1))+1)
cdef complex spectral3(double w, double t): return  2*pi* 2.0 * 6597344572538.565 / (pi * 6024096385542.169 * 2.5448096943108576e-14)  if (w==0) else 2*pi*(2.0*6597344572538.565*6024096385542.169 *w /(pi*(w**2+6024096385542.169**2))) * ((1/(exp((w) * 2.5448096943108576e-14)-1))+1)
cdef complex spectral4(double w, double t): return  2*pi* 2.0 * 6597344572538.565 / (pi * 6024096385542.169 * 2.5448096943108576e-14)  if (w==0) else 2*pi*(2.0*6597344572538.565*6024096385542.169 *w /(pi*(w**2+6024096385542.169**2))) * ((1/(exp((w) * 2.5448096943108576e-14)-1))+1)
cdef complex spectral5(double w, double t): return  2*pi* 2.0 * 6597344572538.565 / (pi * 6024096385542.169 * 2.5448096943108576e-14)  if (w==0) else 2*pi*(2.0*6597344572538.565*6024096385542.169 *w /(pi*(w**2+6024096385542.169**2))) * ((1/(exp((w) * 2.5448096943108576e-14)-1))+1)
cdef complex spectral6(double w, double t): return  2*pi* 2.0 * 6597344572538.565 / (pi * 6024096385542.169 * 2.5448096943108576e-14)  if (w==0) else 2*pi*(2.0*6597344572538.565*6024096385542.169 *w /(pi*(w**2+6024096385542.169**2))) * ((1/(exp((w) * 2.5448096943108576e-14)-1))+1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def cy_td_ode_rhs(
        double t,
        complex[::1] vec,
        complex[::1,:] H0,
        complex[::1,:] A0,
        complex[::1,:] A1,
        complex[::1,:] A2,
        complex[::1,:] A3,
        complex[::1,:] A4,
        complex[::1,:] A5,
        complex[::1,:] A6,
        unsigned int nrows):
    
    cdef double complex * out = <complex *>PyDataMem_NEW_ZEROED(nrows**2,sizeof(complex))
     
    cdef complex[::1, :] H = farray_alloc(nrows)
    cdef complex[::1, :] evecs = farray_alloc(nrows)
    cdef double * eigvals = <double *>PyDataMem_NEW_ZEROED(nrows,sizeof(double))
    dense_add_mult(H, H0, 1)
    ZHEEVR(H, eigvals, evecs, nrows)
    PyDataMem_FREE(&H[0,0])
    cdef double complex * eig_vec = vec_to_eigbasis(vec, evecs, nrows)
    diag_liou_mult(eigvals, eig_vec, out, nrows)
    cdef double[:,::1] skew = <double[:nrows,:nrows]><double *>PyDataMem_NEW_ZEROED(nrows**2,sizeof(double))
    cdef double dw_min = skew_and_dwmin(eigvals, skew, nrows)
    br_term_mult(t, A0, evecs, skew, dw_min, spectral0, eig_vec, out, nrows, 1, 0.1, 1e-12)
    br_term_mult(t, A1, evecs, skew, dw_min, spectral1, eig_vec, out, nrows, 1, 0.1, 1e-12)
    br_term_mult(t, A2, evecs, skew, dw_min, spectral2, eig_vec, out, nrows, 1, 0.1, 1e-12)
    br_term_mult(t, A3, evecs, skew, dw_min, spectral3, eig_vec, out, nrows, 1, 0.1, 1e-12)
    br_term_mult(t, A4, evecs, skew, dw_min, spectral4, eig_vec, out, nrows, 1, 0.1, 1e-12)
    br_term_mult(t, A5, evecs, skew, dw_min, spectral5, eig_vec, out, nrows, 1, 0.1, 1e-12)
    br_term_mult(t, A6, evecs, skew, dw_min, spectral6, eig_vec, out, nrows, 1, 0.1, 1e-12)
    

    cdef np.ndarray[complex, ndim=1, mode='c'] arr_out = vec_to_fockbasis(out, evecs, nrows)
    PyDataMem_FREE(&skew[0,0])
    PyDataMem_FREE(&evecs[0,0])
    PyDataMem_FREE(eigvals)
    PyDataMem_FREE(eig_vec)
    PyDataMem_FREE(out)
    return arr_out
