# location of the Boost Python include files and library
 
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib
 
printRank: printRank.o
	/usr/bin/mpic++ -o printRank  printRank.o -L$(BOOST_LIB) -lboost_mpi -lboost_serialization	
printRank.o: printRank.cpp
	/usr/bin/mpic++ -I$(BOOST_INC) -c printRank.cpp
