#-----------------------------------------------------------------------------
# Replace this with the location of LHAPDF on your system (if desired)
LHAPDFLIB   = /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lhapdf/6.1.6-ikhhed2/lib/

HNNLOHOME        = $(PWD)
SOURCEDIR       = $(PWD)/src
VPATH		= $(DIRS)
BIN		= $(HNNLOHOME)/bin
INCPATH  	= $(SOURCEDIR)/Inc
OUTPUT_OPTION	= -o $(HNNLOHOME)/obj/$@

#   NATIVE -- HNNLO internal routines
#   LHAPDF -- Les Houches library
PDFROUTINES = LHAPDF

#FC = ifort 
#FFLAGS 	=   -O3 -I$(INCPATH)

FC = gfortran
FFLAGS 	=  -fno-automatic -fno-f2c -O1 -I$(INCPATH)


DIRS	=	$(HNNLOHOME):\
		$(HNNLOHOME)/obj:\
		$(SOURCEDIR)/User:$(SOURCEDIR)/Matrix:\
		$(SOURCEDIR)/Need:$(SOURCEDIR)/Phase:\
		$(SOURCEDIR)/Integrate:\
               
# -----------------------------------------------------------------------------
# Specify the object files. 


MATRIXFILES = \
hjetfill.o \
gg_h.o \
gg_hg.o \
gg_hg_rescaled.o \
gg_hg_v.o \
gg_hg_z.o \
gg_hg_gvec.o \
gg_hg_gs.o \
gg_hgg.o \
hqqgg.o \
h4g.o \
h4q.o \
gg_hwwgg.o \
gg_hwwg_gsnew.o \
gg_hwwg_v.o \
gg_hwwg_z.o \
gg_hwwg_gvec.o \
qqb_hww_g.o \
qqb_hww_g_rescaled.o \
qqb_hww.o \
qqb_hzz.o \
qqb_hzz_g.o \
qqb_hzz_g_rescaled.o \
gg_hzzgg.o \
gg_hzzg_gsnew.o \
gg_hzzg_gvec.o \
gg_hzzg_v.o \
gg_hzzg_z.o \
gg_hexact.o \
qqb_hwwexact.o \
qqb_hzzexact.o \


INTEGRATEFILES = \
vegas.o \
mbook.o \
ran0.o \
ran1.o \



NEEDFILES = \
boost.o \
branch.o \
ckmfill.o \
coupling.o \
couplz.o \
dipoles.o \
dipoles_fac.o \
dipolesub.o \
dot.o \
dotem.o \
getptildejet.o \
gtperp.o \
higgsp.o \
higgsw.o \
histofinLH.o \
includedipole.o \
lowintHsp.o \
lowint_incldip.o \
masscuts.o \
hnnlo.o \
hinit.o \
hexit.o \
integrate.o \
ptyrap.o \
r.o \
setup.o \
realint.o \
realvirt2.o \
realvirt4.o \
countH.o \
myli2.o \
myli3.o \
sethparams.o \
smalls.o \
spinork.o \
spinoru.o \
storedip.o \
storeptilde.o \
strcat.o \
swapjet.o \
transform.o \
virtint.o \
writeinfo.o \
zeromsq.o \
besselkfast.o \
isolation.o \
isolation2.o \
isolation4.o \
ddilog.o \
alfamz.o \
newton1.o \
abisq.o \
abisqnnlo.o \
ehsv.o \
li2.o \
c1higgs.o \

PHASEFILES = \
gen2.o \
gen3.o \
gen4MIO.o \
gen4h.o \
gen5.o \
phase5.o \
genBORN2.o \
genBORN4.o \
phase3.o \
phase4.o \
phi1_2.o \
phi3m0.o \
breitw.o \
gen6hp.o \
phase6h.o \
phi1_2h.o \
phi1_2m_bw.o \
phi1_2m.o \



USERFILES = \
deltarj.o \
genclust2.o \
genclust_kt.o \
genclust_cone.o \
cuts.o \
getet.o \
mdata.o \
miscclust.o \
plotter.o \
boostx.o \




LIBDIR=.
LIBFLAGS=

ifeq ($(PDFROUTINES),LHAPDF)
   PDFFILES += \
   pdf_lhapdf.o \
   pdfset_lhapdf.o
   LIBDIR += -L$(LHAPDFLIB)
   LIBFLAGS += -lLHAPDF
   PDFMSG='   ----> HNNLO compiled with LHAPDF routines <----'
else
ifeq ($(PDFROUTINES),NATIVE)
   PDFFILES += \
   pdf.o \
   pdfset.o
   PDFMSG='   ----> HNNLO compiled with its own PDFs <----'
endif
endif

CODE = $(NEEDFILES) $(PHASEFILES) \
          $(USERFILES) $(MATRIXFILES) $(INTEGRATEFILES) $(PDFFILES) \
          

hnnlo: $(CODE)
	$(FC) $(FFLAGS) -L$(LIBDIR) -o $@ \
	$(patsubst %,obj/%,$(CODE)) $(LIBFLAGS)
	mv hnnlo bin/
	@echo $(PDFMSG)

# -----------------------------------------------------------------------------


