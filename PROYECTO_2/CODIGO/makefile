#compilador
FC = gfortran

#directorio subrutinas
SUBDIR=./subrutinas

#SUBRUTINAS = CALCULO_RHS.f90 THOMAS.f90 CC.f90 PREDICCION_VELOCIDAD.f90 \
#CORRECCION_PRESION.f90 CORRECCION_VELOCIDAD.f90 \

all: T2.f90
	$(FC) $(SUBDIR)/CALCULO_RHS.o $(SUBDIR)/THOMAS.o \
	$(SUBDIR)/CC.o $(SUBDIR)/PREDICCION_VELOCIDAD.o \
	$(SUBDIR)/CORRECCION_PRESION.o \
	$(SUBDIR)/CORRECCION_VELOCIDAD.o \
	T2.f90 #-lumfpack -lamd -llapack -lrefblas

compilar_subrutinas:
	cd ./subrutinas && $(MAKE)
