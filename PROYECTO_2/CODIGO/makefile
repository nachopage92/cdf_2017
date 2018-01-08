#compilador
FC = gfortran

#directorio subrutinas
SUBDIR=./subrutinas

#subrutinas
#subrutinas = CALCULO_RHS.o WPENS_U.o WPENS_V.o THOMAS.o CC.o  \
# PREDICCION_VELOCIDAD.o CORRECCION_PRESION.f90 CORRECCION_VELOCIDAD.f90 \
# CALCULO_CONVERGENCIA.o


all: T2.f90
	$(FC) $(SUBDIR)/CALCULO_RHS.o $(SUBDIR)/WPENS_U.o $(SUBDIR)/WPENS_V.o \
	$(SUBDIR)/THOMAS.o $(SUBDIR)/CC.o $(SUBDIR)/PREDICCION_VELOCIDAD.o \
	$(SUBDIR)/TDMA_PHI.o $(SUBDIR)/CORRECCION_PRESION.o \
	$(SUBDIR)/CORRECCION_VELOCIDAD.o $(SUBDIR)/CALCULO_CONVERGENCIA.o \
	$(SUBDIR)/PHI_DIRECTO.o \
	T2.f90 #-lumfpack -lamd -llapack -lrefblas

compilar_subrutinas:
	cd ./subrutinas && $(MAKE)
