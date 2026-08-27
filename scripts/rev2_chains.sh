#!/usr/bin/env bash
# Revisao do paper 2: refazer as cadeias que carregam as afirmacoes fortes,
# TODAS com o perfil `perturbative-gr` (init largo), para que o tratamento seja
# genuinamente identico entre EoS. As seis nao-SLy do paper usaram
# `perturbative-gr-narrowinit` (init_width alpha_U = 0.10), o que nao esta
# declarado no texto e e candidato a sub-exploracao.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
STEPS=1600; WALKERS=32; ENS=8; THREADS=4
run () {  # nome  caminho_eos
  local nm="$1" eos="$2"
  local out="output/rev2/${nm}"
  mkdir -p "$out"
  ./build/brana_tov mcmc $STEPS $WALKERS $ENS $THREADS perturbative-gr "$eos" "$out" \
     > "$out/run.log" 2>&1
  echo "[$(date +%H:%M:%S)] concluida: $nm" >> output/rev2/progress.log
}
run QHC19A  data/eos/tables/QHC19A.txt &
run NL3     data/eos/tables/NL3.txt &
run SLy     data/SLy.txt &
wait
echo "[$(date +%H:%M:%S)] TODAS CONCLUIDAS" >> output/rev2/progress.log
