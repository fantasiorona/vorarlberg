#! /bin/bash
#
g++ -c -Wall ./src/simple_ga.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ simple_ga.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm simple_ga.o
mv a.out SimpleGA
#
echo "Normal end of execution."
