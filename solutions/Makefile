FLAGS = -O3 -DNDEBUG -std=c++11
INCLUDES = /usr/include/eigen3

all: solution1 solution2 solution3 solution4

solution4: solution4.cpp OdeIntTemplate.o
	g++ ${FLAGS} -I ${INCLUDES} -o $@ $^

solution3: solution3.cpp OdeInt.o
	g++ ${FLAGS} -I ${INCLUDES} -o $@ $^

solution2: solution2.cpp
	g++ ${FLAGS} -I ${INCLUDES} -o $@ $^

solution1: solution1.cpp
	g++ ${FLAGS} -I ${INCLUDES} -o $@ $^

%.o: %.cpp
	g++ $(FLAGS) -I ${INCLUDES} -c $^ -o $@ 

clean:
	rm -f *.o solution1 solution2 solution3 solution4

