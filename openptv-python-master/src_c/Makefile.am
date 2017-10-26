## Process this file with automake to produce Makefile.in

check_LTLIBRARIES = libptv.la

libptv_la_SOURCES = tools.c ttools.c imgcoord.c multimed.c \
    ptv.h track.c pointpos.c intersect.c orientation.c epi.c \
    jw_ptv.c segmentation.c peakfitting.c \
    correspondences.c sortgrid.c image_processing.c \
    tracking_run.c tracking_run.h corresp.c btree.c
libptv_la_LIBADD = -lm -loptv 

check_PROGRAMS = main
main_SOURCES = main.c
main_LDADD = libptv.la

