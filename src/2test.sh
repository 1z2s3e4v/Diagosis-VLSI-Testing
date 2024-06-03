#!/bin/bash   
pattern="$1"
circuit="$2"
log="$3"

./../bin_reference/atpg_reference -diag  $pattern $circuit $log  && ./atpg -diag $pattern $circuit $log 