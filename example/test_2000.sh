#!/bin/bash
if [[ $# < 1 ]]
then
	echo Usage: ./test_2000.sh [err]
	exit
fi

err=$1
./szx -z -f -i ~/Data/ARAMCO/new/pressure_2000_1008x1008x352 -3 352 1008 1008 -M REL -R $err
./szx -x -f -i ~/Data/ARAMCO/new/pressure_2000_1008x1008x352 -3 352 1008 1008 -s ~/Data/ARAMCO/new/pressure_2000_1008x1008x352.szx -a
