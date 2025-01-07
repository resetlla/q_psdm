#!/bin/bash

cat ./traveltime/L350.ttt_* > ./traveltime/L350.ttt
cat ./vel/vel_loadin_* > ./vel/vel_loadin

../src/q-psdm path=./0/qmig3dd parfilepath=./0/qmig3dd/par/par.dat ttpath=./traveltime nodename=node2_0