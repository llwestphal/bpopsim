###
##
## AUTHORS
##
## Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
## David B. Knoester
##
## LICENSE AND COPYRIGHT
##
## Copyright (c) 2010 Michigan State University
##
## breseq is free software; you can redistribute it and/or modify it under the terms the 
## GNU General Public License as published by the Free Software Foundation; either 
## version 1, or (at your option) any later version.
##
###

#AM_CPPFLAGS = -I$(top_srcdir)/$(LIBBAM) $(BOOST_CPPFLAGS)
#AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) -static
#LDADD = $(BOOST_PROGRAM_OPTIONS_LIBS) \
#		$(top_srcdir)/$(LIBBAM)/libbam.la \
#		libbreseq.la

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) -static -lgsl
LDADD = $(BOOST_PROGRAM_OPTIONS_LIBS)

bin_PROGRAMS = bpopsim

bpopsim_SOURCES = \
	cPopulation.cc \
	cSubpopulation.cc \
	icsilog.c \
	main.cpp	