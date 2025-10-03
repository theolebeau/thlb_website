#Fortran Compiler

#MPIF= mpif90
HOSTNAME := $(shell hostname)

ifeq ($(HOSTNAME), glx-calcul3)
    F90 = /data/glx-calcul3/data1/tlebeau/anaconda3/bin/mpifort
else ifeq ($(HOSTNAME), glx-localization)
    F90 = /data/glx-localization/tlebeau/anaconda3/bin/mpifort
else
    F90 = /usr/bin/mpifort  # Default MPI compiler
endif

#F90 = /usr/bin/mpifort
#F90_MPI = mpifort
FFLAGS = -O3 -cpp -ffree-line-length-none -std=legacy -mcmodel=large -fno-range-check -g3 #-fstack-protector-all -finit-real=snan -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check
#FFLAGS_MPI = #-I/$home/data/glx-calcul3/data1/tlebeau/openmpi/include -L/$home/data/glx-calcul3/data1/tlebeau/openmpi/lib
#-O3 -cpp -ffree-line-length-none -std=legacy -mcmodel=large -fno-range-check -g3
#-I/$home/data/glx-calcul3/data1/tlebeau/openmpi/include -L/$home/data/glx-calcul3/data1/tlebeau/openmpi/lib
############compileur rdramses#############

rdr: rdramses

rdramses.o: rdramses.f90
	$(F90) $(FFLAGS) -c rdramses.f90

rdramses: rdramses.o
	$(F90) $(FFLAGS) rdramses.o -o rdramses

crdr:
	\rm -f -r rdramses *.o

##########compileur ascii_to_bin ##########

#all: ascii_to_bin

#modules.o: modules.f90
#	$(F90) $(FFLAGS) -c modules.f90

#ascii_to_bin.o:ascii_to_bin.f90
#	$(F90) $(FFLAGS) -c ascii_to_bin.f90

#ascii_to_bin: modules.o ascii_to_bin.o 
#	$(F90) $(FFLAGS) modules.o ascii_to_bin.o -o ascii_to_bin

#clean:
#	\rm -f -r ascii_to_bin *.o

######### compileur amr2map #################

#all: amr2map

#utils.o: utils.f90
#	$(F90) $(FFLAGS) -c utils.f90

#amr2map.o: amr2map.f90
#	$(F90) $(FFLAGS) -c amr2map.f90

#amr2map: utils.o amr2map.o
#	$(F90) $(FFLAGS) utils.o amr2map.o -o amr2map

#clean:
#	\rm -f -r amr2map *.o


##### compileur part2map ######################

#all: part2map

#utils.o: utils.f90
#	$(F90) $(FFLAGS) -c utils.f90

#part2map.o: part2map.f90
#	$(F90) $(FFLAGS) -c part2map.f90

#part2map: utils.o part2map.o
#	$(F90) $(FFLAGS) utils.o part2map.o -o part2map

#clean:
#	\rm -f -r part2map *.o


##########compileur remove_gal ##########

#all: remove_gal

#modules.o: modules.f90
#	$(F90) $(FFLAGS) -c modules.f90

#remove_gal.o:remove_gal.f90
#	$(F90) $(FFLAGS) -c remove_gal.f90

#remove_gal: modules.o remove_gal.o 
#	$(F90) $(FFLAGS) modules.o remove_gal.o -o remove_gal

#clean:
#	\rm -f -r remove_gal *.o


########## main program with all functions available #############

main: mainf90

#open_file.o:open_file.f90
#	$(F90) $(FFLAGS) -c open_file.f90

modules.o: modules.f90
	$(F90) $(FFLAGS) -c modules.f90

mainf90.o: mainf90.f90
	$(F90) $(FFLAGS) -c mainf90.f90

mainf90: modules.o mainf90.o #open_file.o
	$(F90) $(FFLAGS) modules.o mainf90.o -o mainf90 #open_file.o

cmain:
	\rm -f -r mainf90 *.o

#### main with MPI ####

main_mpi: mainf90

#open_file.o:open_file.f90
#	$(F90) $(FFLAGS) -c open_file.f90

modules_mpi.o: modules.f90
	$(F90_MPI) $(FFLAGS_MPI) -c modules.f90

mainf90_mpi.o: mainf90.f90
	$(F90_MPI) $(FFLAGS_MPI) -c mainf90.f90

mainf90_mpi: modules_mpi.o mainf90_mpi.o #open_file.o
	$(F90_MPI) $(FFLAGS_MPI) modules_mpi.o mainf90_mpi.o -o mainf90_mpi #open_file.o

cmain_mpi:
	\rm -f -r mainf90 *.o
