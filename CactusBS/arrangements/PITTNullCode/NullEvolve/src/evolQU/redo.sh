PROG=evolQU
gfortran -o ${PROG} ${PROG}.f90 || exit
rm -f fort.*
./${PROG}
rm -f ./${PROG}
  Xmgrblock  -c 1:3 fort.3*
  Xmgrblock  -c 1:3 fort.1* 


