#compilador
FC = gfortran

#directorio subrutinas
SUBDIR=./subrutinas

all: T2.f90
	$(FC) \
	$(SUBDIR)/VARIABLES.o \
	$(SUBDIR)/CALCULO_RHS.o \
	$(SUBDIR)/THOMAS.o \
	$(SUBDIR)/CONDICIONES_INICIALES.o \
	$(SUBDIR)/CONDICIONES_CONTORNO.o \
	$(SUBDIR)/PREDICCION_VELOCIDAD.o \
	$(SUBDIR)/CORRECCION_PRESION.o \
	$(SUBDIR)/CORRECCION_VELOCIDAD.o \
	$(SUBDIR)/EXPORTAR_DATOS.o \
	$(SUBDIR)/EXPORTAR_DATOS_2.o \
	T2.f90 -o PROYECTO2_IGNACIO_APABLAZA

compilar_subrutinas:
	cd ./subrutinas && $(MAKE)
	
compilar_modulo:
	$(FC) -c VARIABLES.f90
