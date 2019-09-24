
CC = cc
CFLAGS = -O3 -m64 
LDFLAGS= -lm -m64 

cak_objects = cak.o

cak: $(cak_objects)
	$(CC) $(CFLAGS) $(cak_objects) $(LDFLAGS) -o cak
	
clean :
	rm -f *.o  *~

