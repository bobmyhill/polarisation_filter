FC = gfortran
FFLAGS = 
TARGET = polarisation_filter

$(TARGET): four1.o
	${FC} ${FFLAGS} nsyt-filter.f four1.o

four1.o: four1.f
