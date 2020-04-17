CC = mpicc
CFLAGS = -O3 -m64   -I$(PYTHON)/include $(MPI_FLAG)
LDFLAGS= -L$(PYTHON)/lib -lm -lgsl -lgslcblas -m64  

cak_objects = cak.o

MPI_FLAG = -DMPI_ON

cak: $(cak_objects)
	$(CC) $(CFLAGS) $(cak_objects) $(LDFLAGS) -o cak
	
clean :
	rm -f *.o  *~

