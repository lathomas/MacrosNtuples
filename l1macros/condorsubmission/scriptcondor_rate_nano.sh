#!/bin/bash

cd /user/lathomas/GITStuff/L1Nano/new/CMSSW_14_0_1/src
cmsenv
cd /user/lathomas/GITStuff/MacrosNtuples/l1macros

python3 rate_singleobject_nano.py --max_events -1 -i $1 -o $2 


