OBJS=sYlm.o Legendre.o h5read.o
PROGS=h5read

HDF5_DIR=/home/yosef/software/hdf5-1.8.1/
HDF5_LIBS= -lhdf5_hl -lhdf5 -lz
${PROGS}: ${OBJS}
	g++ ${OBJS} -g -o $@  -L ${HDF5_DIR}/lib ${HDF5_LIBS}

%.o: %.cc
	g++ -Wall -O3 -g $< -c -I ${HDF5_DIR}/include

clean:
	rm ${OBJS} ${PROGS}
