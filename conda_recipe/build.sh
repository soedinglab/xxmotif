#!/bin/bash

mkdir build
cd build
cmake .. && make

mkdir -p $PREFIX/bin
cp gimmemotif $PREFIX/bin
ln -s $PREFIX/bin/gimmemotif $PREFIX/bin/XXmotif
