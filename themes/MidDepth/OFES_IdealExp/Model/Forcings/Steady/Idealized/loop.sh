#!/bin/bash

#wind_period="180"
#wind_period="90"
wind_period="140"

for year in `seq 1 30`; do
  echo "year = "${year}
  julia generate.jl "$year" "$wind_period" > log_Year${year}_${wind_period}dy
done

