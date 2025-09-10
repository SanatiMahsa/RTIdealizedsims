#############################################################################
# Compilation time parameters
NVECTOR=32
NDIM=3
NPRE=8
NGROUPS=3
NIONS=3
NCR=0
NENER=0+$(NCR)
NVAR=10+$(NIONS)+$(NENER)
SOLVER=mhd
PATCH= ../patch/IDRTSink_DB
EXEC=ramses_rtmhd-fullphys
ATON_FLAGS= #-DATON  # Uncomment to enable ATON.
#TRACKER_FLAGS= -DTRACKER
DICE_FLAGS= -DDICE
RT_FLAGS= -DRT -DNIONS=$(NIONS) -DNGROUPS=$(NGROUPS)
#############################################################################
# Select line with commons file to use mechanical feedback
MECHFB=mechanical_commons.o
MECHFINE=mechanical_fine.o 
#MECHFB=
#MECHFINE=
# Select line with sink files to use sink particles
SINKFILES=init_sink.o output_sink.o
# SINKFILES=
#############################################################################
# Use write files
WRITEFILES=write_makefile.o write_patch.o write_gitinfo.o 
WRITEFILES=
#############################################################################
# Use conduction and CR files
CONDFILES=conduction.o cond_split.o cond_split_semi.o conduction_fine.o coupling_fine.o
CRFILES=crdiffusion.o crdiff_split.o crdiffusion_fine.o 
CONDFILES=
CRFILES=
#############################################################################
GITBRANCH=$(shell git rev-parse --abbrev-ref HEAD)
GITHASH=$(shell git log --pretty=format:'%H' -n 1)
GITREPO=$(shell git config --get remote.origin.url)
BUILDDATE=$(shell date +"%D-%T")
DEFINES= -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) \
	 -DNENER=$(NENER) -DNCR=$(NCR) -DNVAR=$(NVAR) -DSOLVER$(SOLVER) \
	 $(ATON_FLAGS) $(RT_FLAGS) $(TRACKER_FLAGS) $(DICE_FLAGS)
#############################################################################
# Fortran compiler options and directives
# --- MPI, ifort syntax ------------------------------
F90 = mpif90
F90 = mpiifort
#FFLAGS = -cpp -O3 -xCORE-AVX512 -traceback -DQUADHILBERT $(DEFINES) # -DNOSYSTEM
#FFLAGS = -cpp -O3 -DQUADHILBERT -DLONGINT $(DEFINES) -DNOSYSTEM
#FFLAGS = -cpp -O3 -DQUADHILBERT $(DEFINES) -DNOSYSTEM
FFLAGS = -cpp -O2 -DQUADHILBERT $(DEFINES) -DNOSYSTEM -traceback 
#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI = 
#LIBMPI = -lfmpi -lmpi -lelan
# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart
LIBS = $(LIBMPI)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../$(SOLVER):../aton:../rt:../hydro:../pm:../poisson:../amr
#VPATH = $(PATCH):../cr:../conduction:../$(SOLVER):../aton:../rt:../hydro:../pm:../poisson:../amr
#############################################################################
# All objects
MODOBJ = amr_parameters.o amr_commons.o random.o pm_parameters.o pm_commons.o poisson_parameters.o \
        poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o \
	sparse_mat.o clfind_commons.o gadgetreadfile.o $(MECHFB) $(WRITEFILES) #tracker_commons.o track_center_mass.o 
AMROBJ = read_params.o init_amr.o init_time.o init_refine.o adaptive_loop.o amr_step.o \
	update_time.o output_amr.o flag_utils.o physical_boundaries.o virtual_boundaries.o \
	refine_utils.o nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o \
	units.o light_cone.o #movie.o
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o newdt_fine.o \
	particle_tree.o add_list.o remove_list.o star_formation.o sink_particle.o feedback.o \
	clump_finder.o clump_merger.o flag_formation_sites.o \
	$(MECHFINE) $(SINKFILES)
# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o multigrid_coarse.o \
	multigrid_fine_commons.o multigrid_fine_fine.o multigrid_fine_coarse.o gravana.o \
	boundary_potential.o rho_ana.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.o \
	godunov_fine.o uplmde.o umuscl.o interpol_hydro.o godunov_utils.o condinit.o \
	hydro_flag.o hydro_boundary.o boundana.o read_hydro_params.o synchro_hydro_fine.o \
	$(CONDFILES) $(CRFILES)
# RT objects
RTOBJ = rt_init_hydro.o rt_init_xion.o rt_init.o rt_init_flow_fine.o rt_output_hydro.o \
	rt_godunov_fine.o rt_interpol_hydro.o rt_godunov_utils.o rt_condinit.o \
	rt_hydro_flag.o rt_hydro_boundary.o rt_boundana.o rt_read_hydro_params.o rt_units.o
# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ)
# RT additions to AMR and MOD objects
ifneq (,$(findstring DRT,$(DEFINES)))
MODOBJ += rt_parameters.o rt_hydro_commons.o coolrates_module.o rt_spectra.o rt_cooling_module.o rt_flux_module.o 
AMRLIB += $(RTOBJ)
endif
# ATON objects
ATON_MODOBJ = timing.o radiation_commons.o rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o \
	rad_backup.o ../aton/atonlib/libaton.a
#############################################################################
ramses:	$(MODOBJ) $(AMRLIB) ramses.o
	$(F90) $(MODOBJ) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
	# rm write_makefile.f90
	# rm write_patch.f90
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
	# rm write_makefile.f90
	# rm write_patch.f90
#############################################################################
write_gitinfo.o: FORCE
	$(F90) $(FFLAGS) -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' -DGITHASH=\'"$(GITHASH)"\' \
	-DGITREPO=\'$(GITREPO)\' -DBUILDDATE=\'"$(BUILDDATE)"\' -c ../amr/write_gitinfo.f90 -o $@
write_makefile.o: FORCE
	../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)
	$(F90) $(FFLAGS) -c write_makefile.f90 -o $@
write_patch.o: FORCE
	../utils/scripts/cr_write_patch.sh $(PATCH)
	$(F90) $(FFLAGS) -c write_patch.f90 -o $@
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
FORCE:
#############################################################################
clean :
	rm *.o *.$(MOD)
#############################################################################
