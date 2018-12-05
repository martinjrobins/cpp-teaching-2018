# location of the Boost Python include files and library
 
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib
 
common_factor: common_factor.o
	g++ -o common_factor common_factor.o	
common_factor.o: common_factor.cpp
	g++ -I$(BOOST_INC) -c common_factor.cpp
