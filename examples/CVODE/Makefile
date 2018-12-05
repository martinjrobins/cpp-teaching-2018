CVODE_INC = /usr/include/sundials
CVODE_LIB = /usr/lib
     
cvRoberts_dns: cvRoberts_dns.o
	g++ -O -o $@  $^ -L$(CVODE_LIB) -lsundials_cvode -lsundials_nvecserial 

%.o: %.cpp
	g++ -I$(CVODE_INC) -c $^ -o $@ 

clean:
	rm *.o
	rm cvRoberts_dns
