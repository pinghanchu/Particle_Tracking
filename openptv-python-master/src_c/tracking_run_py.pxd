
cdef extern from "optv/parameters.h":
    ctypedef struct sequence_par:
        int first, last

cdef extern from "tracking_run.h":
    ctypedef struct tracking_run:
        sequence_par *seq_par
    cdef void tr_free(tracking_run *tr)

cdef class TrackingRun:
    cdef tracking_run *tr

