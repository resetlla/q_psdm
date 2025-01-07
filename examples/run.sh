#!/bin/bash

cat ./traveltime/L350.ttt_* > ./traveltime/L350.ttt
cat ./vel/vel_loadin_* > ./vel_vel_loadin

../src/qmig3dd_multithreading path=./0/qmig3dd parfilepath=./0/qmig3dd/par/par.dat ttpath=./traveltime nodename=node2_0