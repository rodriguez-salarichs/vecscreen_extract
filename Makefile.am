man_MANS = vecscreen_statistic.man
bin_PROGRAMS = vecscreen_statistic 
vecscreen_statistic_SOURCES = vecscreen_statistic.f
AM_vecscreen_statistic_OBJECTS = vecscreen_statistic.o
AM_FFLAGS=-w

vecscreen_statistic.o: vecscreen_statistic.f
	$(F77) -w -c vecscreen_statistic.f

EXEEXT = .x
