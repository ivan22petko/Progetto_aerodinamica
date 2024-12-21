#!/bin/bash

chmod +x run_xfoil_2.sh

IFILE='instruction_xfoil.txt'

echo "LOAD BL207.dat" > $IFILE
echo "OPER" >> $IFILE
echo "VISC" >> $IFILE
echo "25000000" >> $IFILE
echo "ALFA 2" >> $IFILE
echo "VPLO" >> $IFILE
echo "CF" >> $IFILE
echo "DUMP dati_6.dat" >> $IFILE
echo "plCf" >> $IFILE
echo "HARD cfplot.png" >> $IFILE
echo '  ' >> $IFILE
echo "QUIT" >> $IFILE
xfoil < $IFILE
