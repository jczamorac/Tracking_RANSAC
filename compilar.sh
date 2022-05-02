#!/bin/sh

cd build
##rm -r *
cmake ../
make
cp libTrackinglib.so ../test/lib/
cd ../include
cp *.h ../test/include/
cd ../
