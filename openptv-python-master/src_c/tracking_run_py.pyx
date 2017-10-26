# Provides a Python wrapper to the tracking_run struct. This should allow a
# tracking algorithm to be managed from Python, with only the bottleneck parts
# running in C. The tracking_run struct holds information needed or generated
# by the different tracking components.
    
from libc.stdlib cimport free

cdef class TrackingRun:
    def get_sequence_range(TrackingRun self):
        return self.tr[0].seq_par[0].first, self.tr[0].seq_par[0].last

    def __dealloc__(TrackingRun self):
        tr_free(self.tr)
        free(self.tr)

