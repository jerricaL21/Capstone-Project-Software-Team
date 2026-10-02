#!/bin/bash
# Davis Lab — staged NAMD launcher (GPU). Run stages IN ORDER.
set -e
NAMD=namd3   # or your cluster's namd binary/module
$NAMD +p8 +devices 0 01_minimize.namd  > 01_minimize.log
$NAMD +p8 +devices 0 02_heat.namd      > 02_heat.log
$NAMD +p8 +devices 0 03_equil.namd     > 03_equil.log
$NAMD +p8 +devices 0 04_production.namd > 04_production.log
echo "Done. Run: python analyze_distance.py"
