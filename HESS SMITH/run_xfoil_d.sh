#!/bin/bash

chmod +x run_xfoil_2.sh

IFILE='instruction_xfoil.txt'

echo "NACA 0012" > $IFILE
echo "PPAR" >> $IFILE
echo "N 102" >> $IFILE
echo "" >> $IFILE
echo "" >> $IFILE  
echo "OPER" >> $IFILE
echo "ALFA 2" >> $IFILE
echo "CPWR naca_0012_cp.txt" >> $IFILE 
echo '  ' >> $IFILE
echo "QUIT" >> $IFILE
xfoil < $IFILE
echo '' 
 
 
 
echo "NACA 0012" > $IFILE
echo "PPAR" >> $IFILE
echo "N 102" >> $IFILE
echo "" >> $IFILE
echo "" >> $IFILE  
echo "OPER" >> $IFILE
echo "PACC" >> $IFILE
echo "polar_cl_d.txt" >> $IFILE
echo "" >> $IFILE 
echo "ASEQ" >> $IFILE
echo "0" >> $IFILE
echo "5" >> $IFILE
echo "0.5" >> $IFILE
echo '  ' >> $IFILE
echo "QUIT" >> $IFILE
xfoil < $IFILE

echo "NACA 0012" > $IFILE
echo "PPAR" >> $IFILE
echo "N 102" >> $IFILE
echo "" >> $IFILE
echo "" >> $IFILE  
echo "OPER" >> $IFILE
echo "PACC" >> $IFILE
echo "polar_cm_d.txt" >> $IFILE
echo "" >> $IFILE 
echo "ASEQ" >> $IFILE
echo "0" >> $IFILE
echo "2" >> $IFILE
echo "0.25" >> $IFILE
echo '  ' >> $IFILE
echo "QUIT" >> $IFILE
xfoil < $IFILE






    
