#! /bin/bash

make -j;
for n in 1000;
do
./main $n 0 1;
done

