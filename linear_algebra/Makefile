all: use_vectors

#For debugging
OPT=-g
#For optimistaion
#OPT=-O

#All objects (except main) come from cpp and hpp 
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<
#use_vectors relies on objects which rely on headers
use_vectors:	use_vectors.cpp Vector.o Exception.o
		g++ ${OPT} -o use_vectors use_vectors.cpp Vector.o Exception.o

clean:
	rm -f *.o *~ use_vectors
