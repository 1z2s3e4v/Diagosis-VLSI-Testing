#!/bin/bash   
pattern="$1"
circuit="$2"
log="$3"

./atpg -diag $pattern $circuit $log 
