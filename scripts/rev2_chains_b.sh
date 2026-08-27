#!/usr/bin/env bash
# Segunda leva: as quatro EoS restantes com init largo, para que o tratamento seja
# identico as tres ja refeitas. Necessario porque o init estreito sub-explorou
# QHC19-A, e o mesmo pode valer para estas.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
STEPS=1600; WALKERS=32; ENS=8; THREADS=3
run () { local nm="$1" eos="$2"; local out="output/rev2/${nm}"; mkdir -p "$out"
  ./build/brana_tov mcmc $STEPS $WALKERS $ENS $THREADS perturbative-gr "$eos" "$out" > "$out/run.log" 2>&1
  echo "[$(date +%H:%M:%S)] concluida: $nm" >> output/rev2/progress.log; }
run BSK24  data/eos/tables/BSK24.txt &
run DD2    data/eos/tables/DD2.txt &
run QHC21D data/eos/tables/QHC21D.txt &
run SPGM4  data/eos/tables/SPGM4.txt &
wait
echo "[$(date +%H:%M:%S)] SEGUNDA LEVA CONCLUIDA" >> output/rev2/progress.log
