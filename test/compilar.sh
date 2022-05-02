#!/bin/sh

cd build
##rm -r *
cmake ../
make
mv find_tracks ../
cd ../
