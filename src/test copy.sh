#!/bin/bash   
pattern="$1"
circuit="$2"
wire="$3"
gate="$4"
inout="$5"
fault="$6"
echo "Reference binary output: "
../bin_reference/atpg_reference  -genFailLog $pattern $circuit -fault $wire $gate $inout $fault
echo "My output: "
./atpg -genFailLog $pattern $circuit -fault $wire $gate $inout $fault