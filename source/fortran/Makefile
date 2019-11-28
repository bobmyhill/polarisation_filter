FC = gfortran
FFLAGS =
TARGET = polarisation_filter

$(TARGET): four1.o
	${FC} ${FFLAGS} filter.f four1.o /usr/local/sac/lib/sacio.a -o $(TARGET)

four1.o: four1.f
