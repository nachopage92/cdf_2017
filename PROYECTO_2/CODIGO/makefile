#compilador
FC = gfortran

#directorio subrutinas
SUBDIR=./subrutinas

#subrutinas
#subrutinas = CALCULO_RHS.o WPENS.o THOMAS.o CC.o RIGHT_HAND_SIDE.o PREDICCION_VELOCIDAD.o


all: T2.f90
	$(FC) $(SUBDIR)/CALCULO_RHS.o $(SUBDIR)/WPENS.o $(SUBDIR)/THOMAS.o	\
	$(SUBDIR)/CC.o $(SUBDIR)/RIGHT_HAND_SIDE.o $(SUBDIR)/PREDICCION_VELOCIDAD.o \
	$(SUBDIR)/TDMA_PHI.o \
	 T2.f90

compilar_subrutinas:
	cd ./subrutinas && $(MAKE)
