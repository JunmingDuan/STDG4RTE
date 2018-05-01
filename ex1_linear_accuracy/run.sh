#! /bin/bash

make -j;
#for n in 4 8 16 24 32 40 48 56 64 72 80 512 1024 2048 4096;
#for n in 10 20 40 80 160 320 640 1280; #  2560 5120;
for n in 10 20 40 80 160 320 640;
#for n in 4;
do
./main $n 0 1;
done

