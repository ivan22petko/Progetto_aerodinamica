#!/bin/bash

chmod +x run_xfoil_2.sh

IFILE='instruction_xfoil.txt'

echo "LOAD BL207.dat" > $IFILE
echo "OPER" >> $IFILE
echo "ALFA 1.276391" >> $IFILE
echo "CPWR cp_th.txt" >> $IFILE
echo "plCp" >> $IFILE
echo "HARD cp_plot.png" >> $IFILE
echo '  ' >> $IFILE
echo "QUIT" >> $IFILE
xfoil < $IFILE

