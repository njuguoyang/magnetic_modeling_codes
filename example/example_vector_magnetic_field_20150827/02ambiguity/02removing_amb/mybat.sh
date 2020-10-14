#!/bin/sh

for file in $(ls field*.dat); do
  echo $file > par1
  awk 'match($0,"field") == 0 {print $0}' par >> par1
  mv par1 par
  ./ambig
  str1=`echo $file | awk -F[d] '{print $2}'`
  str2=$str1'dat'
  mv azimuth.dat 'azimuth'$str2
done
exit 0
