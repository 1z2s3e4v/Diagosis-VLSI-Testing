#!/bin/bash   

./atpg -genFailLog ../patterns/golden_c17.ptn ../sample_circuits/c17.ckt -fault 23GAT"("9")" g5 GO SA1 > log
../bin_reference/atpg_reference -genFailLog ../patterns/golden_c17.ptn ../sample_circuits/c17.ckt -fault 23GAT"("9")" g5 GO SA1 > log2
diff log log2
