# location of the Boost include files and library
 
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib
 
harmonic_oscillator: harmonic_oscillator.o
	g++ -o $@ -L$(BOOST_LIB) $^

%.o: %.cpp
	g++ -I$(BOOST_INC) -c $^ -o $@ 

clean:
	rm *.o
	rm harmonic_oscillator
