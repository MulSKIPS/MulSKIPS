# Macro-Definitions:

#what to build
#DESTDIR=../../../bin
DESTDIR = .
MODDIR  = $(DESTDIR)/modules
INOUDIR = $(DESTDIR)/inpout
INITDIR = $(DESTDIR)/initial
MCEVDIR = $(DESTDIR)/mcevent
LATTDIR = $(DESTDIR)/lattice

LIB= -L /opt/intel/mkl/10.0.1.014/lib/em64t -lmkl_solver -lmkl_lapack -lmkl_em64t -lmkl -lguide -lpthread
what: $(DESTDIR)/mulskips.e

SOURCE = $(MODDIR)/defdertype.f $(MODDIR)/defsystem.f $(MODDIR)/definitions.f \
$(INITDIR)/SetAPB.f  $(INITDIR)/SetInvPC.f    $(INITDIR)/SetInvPwAPB.f $(INITDIR)/SetNCSp.f $(INITDIR)/SetNCSpSi.f \
$(INITDIR)/SetAPBSym.f $(INITDIR)/SetInvP.f  $(INITDIR)/SetInvPwOnlyCorSi.f  $(INITDIR)/SetSiC3C.f \
$(INITDIR)/SetSiCSF.f $(INITDIR)/SetSiGex.f \
$(INITDIR)/SetSiC3Cww.f $(INITDIR)/SetSiC3Cwwaperture.f $(INITDIR)/SetSi.f $(INITDIR)/SetSiFINFET.f \
$(INITDIR)/SetCAD.f $(INITDIR)/SetT.f \
$(INITDIR)/SetAPBSymZimb1.f  $(INITDIR)/SetInvPOff.f  $(INITDIR)/SetInvPwOnlySiorC.f \
$(INITDIR)/SetTrench.f $(INITDIR)/SetAPBSymZimb.f   $(INITDIR)/SetInvPSi.f   $(INITDIR)/SetNC.f \
$(INOUDIR)/FileClose.f  $(INOUDIR)/GetOutputFileName.f  $(INOUDIR)/WriteMolMolXYZFile.f $(INOUDIR)/FileOpen.f $(INOUDIR)/WriteMolMolSource.f \
$(LATTDIR)/FIND_NN.f  $(LATTDIR)/FIND_SN.f \
$(MCEVDIR)/Add.f $(MCEVDIR)/Evaporation.f $(MCEVDIR)/PickTreeEvent.f  $(MCEVDIR)/UpdateTree.f \
$(MCEVDIR)/AllocateArrays.f  $(MCEVDIR)/Get_Prob.f  $(MCEVDIR)/random.f \
$(MCEVDIR)/deposition.f  $(MCEVDIR)/Get_Prob_Ini.f $(MCEVDIR)/SetProbability.f\
$(MCEVDIR)/Erase.f  $(MCEVDIR)/IF_Arm_Choice.f $(MCEVDIR)/UpdateCoor.f main.f \
$(MCEVDIR)/absorption.f  $(MCEVDIR)/desorption.f 

OBJS = $(SOURCE:.f=.o)
clean :         
	rm $(OBJS) mulskips.e *.mod
       
LINTFORSRC = $(SOURCE:.f=.ln)
FC=f95
#FC=gfortran

LINT = lintfor

FFLAGS= -O4 -mcmodel=medium  	
#FFLAGS= -O4 -mcmodel=medium -fcheck=bounds 	

#FFLAGS= -O4 -mcmodel medium -shared-intel	
#FFLAGS= -O5 -pg
#FFLAGS= -C -warn argument_checking -O4 -std90 -warn uninitialized
#FFLAGS= -C -warn argument_checking -O4 -std90 -g -check underflow -check overflow
#FFLAGS= -O4 -g -std90 -C -check underflow -check overflow -check format \
# -check output_conversion -warn argument_checking -warn uninitialized
# Commands:

$(DESTDIR)/mulskips.e:	$(OBJS)
	@echo single file...
	@echo Now compiling... 
	$(FC)  $(FFLAGS) $(OBJS) -o $(DESTDIR)/mulskips.e # $(LIB) 
#	$(FC)  $(FFLAGS) $(OBJS) -o $(DESTDIR)/code14.e $(LIB) -lm -lc
#	$(FC)  $(FFLAGS) $(OBJS) -o $(DESTDIR)/code14.e $(LIB) -lm -lc -L/usr/lib/libphpfp.a \
#	-L/usr/lib/libshpf.a -L/usr/lib/libphpf.a -L/usr/lib/libc.a 
