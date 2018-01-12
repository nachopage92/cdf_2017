#compilador
FC = gfortran

#directorio subrutinas
SUBDIR=./subrutinas

all: T2.f90
	$(FC) $(SUBDIR)/CALCULO_RHS.o \
	$(SUBDIR)/THOMAS.o $(SUBDIR)/CC.o $(SUBDIR)/PREDICCION_VELOCIDAD.o \
	$(SUBDIR)/CORRECCION_PRESION.o \
	$(SUBDIR)/CORRECCION_VELOCIDAD.o $(SUBDIR)/CALCULO_CONVERGENCIA.o \
	$(SUBDIR)/PHI_DIRECTO.o \
	T2.f90 #-lumfpack -lamd -llapack -lrefblas

compilar_subrutinas:
	cd ./subrutinas && $(MAKE)
