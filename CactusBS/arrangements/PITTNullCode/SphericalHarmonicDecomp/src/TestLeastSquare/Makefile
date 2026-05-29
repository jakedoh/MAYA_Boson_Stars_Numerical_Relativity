OBJS1=decomp.o matrix.o sYlm.o driver.o Legendre.o Chebyshev.o
OBJS2=matrix.o testmatrix.o
PROGS=decomp testmatrix
GSL_DIR=~/local/gsl

all: ${PROGS}
decomp: ${OBJS1}
	g++ ${OBJS1} -g -o $@ -L ${GSL_DIR}/lib/  -lgsl -lgslcblas -lm

testmatrix: ${OBJS2}
	g++ ${OBJS2} -g -o $@  -L ${GSL_DIR}/lib/  -lgsl -lgslcblas -lm

%.o: %.cc
	g++ -Wall -O3 -g $< -c -I ${GSL_DIR}/include -DDEBUG_

clean:
	-rm ${OBJS1} ${OBJS2} ${PROGS}
