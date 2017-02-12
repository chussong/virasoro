AUTOMAKE_OPTIONS =
bin_PROGRAMS = virasoro
libexec_PROGRAMS = $(VWSTP)
EXTRA_PROGRAMS = vwstp
LDADD = $(AM_LDADD) @staticlibs@
IDIR = include/
virasoro_SOURCES = main.cpp virasoro.cpp runfile.cpp access.cpp $(IDIR)virasoro.h $(IDIR)runfile.h $(IDIR)hmn.h $(IDIR)cpqmn.h $(IDIR)access.h $(IDIR)mpreal.h $(IDIR)mpcomplex.h
virasoro_LDADD = $(LDADD)
vwstp_SOURCES = vwstp.cpp virasoro.cpp runfile.cpp access.cpp $(IDIR)virasoro.h $(IDIR)runfile.h $(IDIR)hmn.h $(IDIR)cpqmn.h $(IDIR)access.h $(IDIR)mpreal.h $(IDIR)mpcomplex.h
nodist_vwstp_SOURCES = vwstptm.c
vwstp_CFLAGS = $(AM_CFLAGS) @wstpflags@
vwstp_CXXFLAGS = $(AM_CXXFLAGS) @wstpflags@
vwstp_CPPFLAGS = -DHAVE_WSTP_
vwstp_LDADD = $(LDADD) @wstpextralibs@
PKGDIR = @pkgdatadir@

if MAKE_VWSTP
pkgdata_DATA = Virasoro.m
BUILT_SOURCES = vwstptm.c
CLEANFILES = vwstptm.c
vwstptm.c: vwstp.tm
	$(WSPREP) $? -o $@

install-data-hook:
	chown $(USERNAME) $(PKGDIR)/Virasoro.m
endif
