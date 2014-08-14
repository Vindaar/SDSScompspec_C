# Makefile for SDSS-spec
NAME = SDSS-spec
PROGNAME = SDSScompspec

# Linux
CC = gcc
FC = gfortran
CFLAGS = -g -O0 -Wall -D__LINUX -DLINUX -I./ -I/usr/lib/ \
         -I/usr/local/pgplot/ \
	 -L/usr/local/pgplot/ -lpgplot
FFLAGS = -O3 \
         -I/usr/local/pgplot/ \
	 -L/usr/local/pgplot/ -lpgplot
LIBS += -lm -lcpgplot -lpgplot -lcfitsio -lX11 -ldl  \
         -I/usr/local/pgplot/ \
	 -L/usr/local/pgplot/
#TARGET = ${HOME}/usr/vincent/SDSS-compspec/${PROGNAME}
TARGET = ./bin/${PROGNAME}
LINK = ./${PROGNAME}
#LINK = ${HOME}/usr/vincent/SDSS-compsec/${PROGNAME}
# Solaris
#CC = gcc
#FC = f77
#CFLAGS = -O2 -Wall -I/home/mim/include -I./
#LIBS = -lm -L/home/mim/lib -lcpgplot -lpgplot -lcfitsio -lX11 -lsocket -lnsl
#TARGET = /home/mim/bin/SDSS-spec

#ast_epoch2jd.o ast_eq2gal.o ast_precess.o ast_rotmatrix.o 

OBJECTS = $(NAME).o djmax.o djmin.o dust_extinct.o faskropen.o		\
faskwopen.o fcompl.o fjmax.o fjmin.o get_input.o getscbc.o idxdval.o	\
idxival.o ijmax.o isdir.o isodd.o median.o medianrun.o pg_get_wins.o	\
pg_open.o pg_win_rename.o pythag.o qsort_dbleint.o qsort_darray.o	\
res_func.o SDSS_CIVbalnicity.o SDSS_Ebv.o SDSS_compfits.o		\
SDSS_compspec.o SDSS_compinit.o SDSS_Gal_extinct_correct.o		\
SDSS_pgenv_init.o SDSS_primtarg.o SDSS_plot_spec.o SDSS_red_powerlaw.o	\
SDSS_rinputfile.o SDSS_rspSpec.o SDSS_rplate.o SDSS_rspec.o SDSS_restinit.o \
SDSS_restspec.o SDSS_statistics.o stats.o svbksb.o svdcmp.o svdfit.o	\
svdfit_poly.o svdvar.o fitserrmsg.o SDSS_colors.o memory.o error.o      \
ast_eq2gal.o ast_precess.o ast_rotmatrix.c ast_epoch2jd.o


DUST-OBJECTS = subs_lambert.o subs_fits.o subs_memory.o subs_inoutput.o

$(NAME): $(OBJECTS) $(DUST-OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(DUST-OBJECTS) $(LIBS)

link:
	ln -s $(TARGET) $(LINK)

depend:
	makedepend -f Makefile -Y -- $(CFLAGS) -- -s "# Dependencies" \
	$(OBJECTS:.o=.c) $(FORT-OBJECTS:.o=.f)# >& /dev/null

clean: 
	rm -f *~ *.o
	rm -f ./bin/*

all:
	make clean
	make
	make link

# Dependencies

SDSS-spec.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS-spec.o: ../include/stats.h ../include/const.h
SDSS-spec.o: ../include/file.h ../include/fit.h
SDSS-spec.o: ../include/memory.h ../include/error.h
djmax.o: ../include/error.h
djmin.o: ../include/error.h
dust_extinct.o: ../include/error.h
dust_extinct.o: ../include/dust_extinct.h
faskropen.o: ../include/file.h input.h
faskropen.o: ../include/error.h
faskwopen.o: ../include/file.h input.h
faskwopen.o: ../include/error.h
fcompl.o: charstr.h ../include/file.h
fcompl.o: ../include/error.h
fjmax.o: ../include/error.h
fjmin.o: ../include/error.h
get_input.o: charstr.h ../include/error.h input.h
getscbc.o: charstr.h input.h ../include/error.h
ijmax.o: ../include/error.h
median.o: ../include/sort.h ../include/stats.h
median.o: ../include/memory.h ../include/error.h
medianrun.o: ../include/stats.h ../include/sort.h
medianrun.o: ../include/memory.h ../include/error.h
pg_get_wins.o: pg_plot.h charstr.h utils.h ../include/error.h
pg_open.o: pg_plot.h charstr.h utils.h ../include/error.h
pg_win_rename.o: pg_plot.h charstr.h utils.h ../include/error.h
pythag.o: ../include/fit.h
qsort_dbleint.o: ../include/sort.h
res_func.o: ../include/error.h
SDSS_CIVbalnicity.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_CIVbalnicity.o: ../include/stats.h
SDSS_CIVbalnicity.o: ../include/memory.h
SDSS_CIVbalnicity.o: ../include/const.h
SDSS_CIVbalnicity.o: ../include/error.h
SDSS_Ebv.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_Ebv.o: ../include/memory.h ../include/error.h
SDSS_Ebv.o: ../include/interface.h
SDSS_Ebv.o: ../include/subs_fits.h
SDSS_Ebv.o: ../include/subs_lambert.h
SDSS_compfits.o: SDSS-spec.h pg_plot.h charstr.h utils.h astron.h
SDSS_compfits.o: ../include/const.h
SDSS_compfits.o: ../include/memory.h
SDSS_compfits.o: ../include/error.h
SDSS_compspec.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_compspec.o: ../include/const.h
SDSS_compspec.o: ../include/memory.h
SDSS_compspec.o: ../include/error.h
SDSS_compinit.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_compinit.o: ../include/const.h
SDSS_compinit.o: ../include/memory.h
SDSS_compinit.o: ../include/error.h
SDSS_Gal_extinct_correct.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_Gal_extinct_correct.o: ../include/error.h
SDSS_pgenv_init.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_primtarg.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_primtarg.o: ../include/error.h
SDSS_plot_spec.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_plot_spec.o: ../include/memory.h
SDSS_plot_spec.o: ../include/error.h
SDSS_red_powerlaw.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_red_powerlaw.o: ../include/fit.h
SDSS_red_powerlaw.o: ../include/stats.h
SDSS_red_powerlaw.o: ../include/const.h
SDSS_red_powerlaw.o: ../include/memory.h
SDSS_red_powerlaw.o: ../include/error.h
SDSS_rinputfile.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_rinputfile.o: ../include/file.h
SDSS_rinputfile.o: ../include/const.h
SDSS_rinputfile.o: ../include/error.h
SDSS_rplate.o: SDSS-spec.h pg_plot.h charstr.h utils.h astron.h
SDSS_rplate.o: ../include/const.h
SDSS_rplate.o: ../include/memory.h
SDSS_rplate.o: ../include/error.h
SDSS_rspSpec.o: SDSS-spec.h pg_plot.h charstr.h utils.h astron.h
SDSS_rspSpec.o: ../include/const.h
SDSS_rspSpec.o: ../include/memory.h
SDSS_rspSpec.o: ../include/error.h
SDSS_rspec.o: SDSS-spec.h pg_plot.h charstr.h utils.h astron.h
SDSS_rspec.o: ../include/const.h
SDSS_rspec.o: ../include/memory.h
SDSS_rspec.o: ../include/error.h
SDSS_restinit.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_restinit.o: ../include/memory.h
SDSS_restinit.o: ../include/error.h
SDSS_restspec.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_restspec.o: ../include/const.h
SDSS_restspec.o: ../include/memory.h
SDSS_restspec.o: ../include/error.h
SDSS_statistics.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_statistics.o: ../include/stats.h
SDSS_statistics.o: ../include/error.h
SDSS_statistics.o: ../include/memory.h
stats.o: ../include/stats.h ../include/error.h
svbksb.o: ../include/memory.h ../include/error.h
svdcmp.o: ../include/fit.h ../include/memory.h
svdcmp.o: ../include/error.h
svdfit.o: ../include/fit.h ../include/memory.h
svdfit.o: ../include/error.h
svdvar.o: ../include/memory.h ../include/error.h
fitserrmsg.o: ../include/usual.h ../include/error.h
fitserrmsg.o: ../include/memory.h ../include/sort.h
SDSS_colors.o: SDSS-spec.h pg_plot.h charstr.h utils.h
SDSS_colors.o: ../include/error.h colorcurves.h
SDSS_colors.o: ../include/memory.h
memory.o: ../include/memory.h ../include/error.h
ast_eq2gal.o: astron.h ../include/const.h
ast_eq2gal.o: ../include/error.h
ast_precess.o: astron.h ../include/const.h
ast_precess.o: ../include/memory.h
ast_precess.o: ../include/error.h
ast_rotmatrix.o: astron.h ../include/const.h
ast_epoch2jd.o: astron.h ../include/const.h
