#!/bin/bash

cd C/library/source
make

cd ../../data_simulation/source
make

cd ../../eta_simulation
./compile

cd ../
