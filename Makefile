SHELL := /bin/bash
c_compiler = gcc#Select the C compiler
###Compiling option###
#GSL
gsl = -lgsl -lgslcblas -lm
gflags = -ffast-math
###Directories###
#Objects#
o_dir = ./objects/
####  Structures   ####
s_dir = ./structures/
#Systems#
#3D#
HS = Hard_Sphere
HSSW = $(HS)_Square_Well
HSDBLEXP = $(HS)_Double_Exp
HSDBLEYUK = $(HS)_Double_Yukawa
#2D#
HD = Hard_Disk
s_libraries = $(HS) $(HSSW) $(HD) $(HSDBLEXP) $(HSDBLEYUK)
#####  Dynamics  #####
d_dir = ./dynamics/
#Dynamic libraries#
SCGLE = SCGLE
NESCGLE = NESCGLE
d_libraries = $(SCGLE) $(NESCGLE)
##### MATH ######
m_dir = ./math/
#Math libraries#
math_aux = math_aux
m_libraries = $(math_aux)
#all_dep = structure_test.c
all_dep = test.c \
$(s_dir)structures.h $(s_dir)structures.c $(o_dir)structures.o\
$(foreach str,$(s_libraries),$(o_dir)$(str).o)\
$(foreach str,$(s_libraries),$(s_dir)$(str)/$(str).c)\
$(foreach str,$(s_libraries),$(s_dir)$(str)/$(str).h)\
$(d_dir)/dynamics.h\
$(foreach str,$(d_libraries),$(o_dir)$(str).o)\
$(foreach str,$(d_libraries),$(d_dir)$(str)/$(str).c)\
$(foreach str,$(d_libraries),$(d_dir)$(str)/$(str).h)\
$(m_dir)math_aux.c $(m_dir)math_aux.h $(o_dir)math_aux.o

test.exe: $(all_dep)
	gcc $(o_dir)*.o -o test.exe test.c $(gsl) $(gflags)


#Building structural libraries
structures_dep = $(s_dir)structures.c $(s_dir)structures.h \
$(foreach str,$(s_libraries),$(o_dir)$(str).o)\
$(foreach str,$(s_libraries),$(s_dir)$(str)/$(str).c)\
$(foreach str,$(s_libraries),$(s_dir)$(str)/$(str).h)
$(o_dir)structures.o: $(structures_dep)
	gcc -c $(s_dir)structures.c -o $(o_dir)structures.o $(gsl) $(gflags)

HSdir = $(s_dir)$(HS)/
HSdep = $(HSdir)$(HS).c $(HSdir)$(HS).h
$(o_dir)$(HS).o:  $(HSdep)
	gcc -c $(HSdir)$(HS).c -o $(o_dir)$(HS).o $(gsl) $(gflags)

HSSWdir = $(s_dir)$(HSSW)/
HSSWdep = $(HSSWdir)$(HSSW).c $(HSSWdir)$(HSSW).h $(o_dir)$(HS).o
$(o_dir)$(HSSW).o: $(HSSWdep)
	gcc -c $(HSSWdir)$(HSSW).c -o $(o_dir)$(HSSW).o $(gsl) $(gflags)

HSDBLEXPdir = $(s_dir)$(HSDBLEXP)/
HSDBLEXPdep = $(HSDBLEXPdir)$(HSDBLEXP).c $(HSDBLEXPdir)$(HSDBLEXP).h $(o_dir)$(HS).o
$(o_dir)$(HSDBLEXP).o: $(HSDBLEXPdep)
	gcc -c $(HSDBLEXPdir)$(HSDBLEXP).c -o $(o_dir)$(HSDBLEXP).o $(gsl) $(gflags)

HSDBLEYUKdir = $(s_dir)$(HSDBLEYUK)/
HSDBLEYUKdep = $(HSDBLEYUKdir)$(HSDBLEYUK).c $(HSDBLEYUKdir)$(HSDBLEYUK).h $(o_dir)$(HS).o
$(o_dir)$(HSDBLEYUK).o: $(HSDBLEYUKdep)
	gcc -c $(HSDBLEYUKdir)$(HSDBLEYUK).c -o $(o_dir)$(HSDBLEYUK).o $(gsl) $(gflags)

HDdir = $(s_dir)$(HD)/
HDdep = $(HDdir)$(HD).c $(HDdir)$(HD).h
$(o_dir)$(HD).o: $(HDdep)
	gcc -c $(HDdir)$(HD).c -o $(o_dir)$(HD).o $(gsl) $(gflags)

#Building dynamical libraries
SCGLEdir = $(d_dir)$(SCGLE)/
SCGLEdep = $(SCGLEdir)$(SCGLE).c $(SCGLEdir)$(SCGLE).h
$(o_dir)$(SCGLE).o: $(SCGLEdep)
	gcc -c $(SCGLEdir)$(SCGLE).c -o $(o_dir)$(SCGLE).o $(gsl) $(gflags)

NESCGLEdir = $(d_dir)$(NESCGLE)/
NESCGLEdep = $(NESCGLEdir)$(NESCGLE).c $(NESCGLEdir)$(NESCGLE).h
$(o_dir)$(NESCGLE).o: $(NESCGLEdep)
	gcc -c $(NESCGLEdir)$(NESCGLE).c -o $(o_dir)$(NESCGLE).o $(gsl) $(gflags)

#Building math libraries
math_auxdir = $(m_dir)
math_auxdep = $(math_auxdir)$(math_aux).c $(math_auxdir)$(math_aux).h
$(o_dir)$(math_aux).o: $(math_auxdep)
	gcc -c $(math_auxdir)$(math_aux).c -o $(o_dir)$(math_aux).o $(gsl) $(gflags)



clean:
	rm $(o_dir)*.o
	rm ./*.dat
	rm ./*.exe
