TESTS = check_ptv check_epi check_corresp 
XFAIL_TESTS = check_corresp
check_PROGRAMS = check_ptv check_epi check_corresp

check_ptv_SOURCES = check_ptv.c $(top_builddir)/src_c/ptv.h
check_ptv_CFLAGS = @CHECK_CFLAGS@ 
check_ptv_LDADD = $(top_builddir)/src_c/libptv.la @CHECK_LIBS@ -lm -loptv

check_epi_SOURCES = check_epi.c $(top_builddir)/src_c/calibration.h
check_epi_CFLAGS = @CHECK_CFLAGS@
check_epi_LDADD = $(top_builddir)/src_c/libptv.la @CHECK_LIBS@ -lm -loptv

check_corresp_SOURCES = check_corresp.c $(top_builddir)/src_c/calibration.h
check_corresp_CFLAGS = @CHECK_CFLAGS@
check_corresp_LDADD = $(top_builddir)/src_c/libptv.la @CHECK_LIBS@ -lm

