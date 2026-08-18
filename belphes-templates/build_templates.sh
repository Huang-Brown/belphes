#!/usr/bin/env bash
#
# build_templates.sh -- driver for steps 1-2 of the belphes template pipeline.
#
# Compiles the tools, optionally fetches NanoAOD from CERN Open Data, builds the
# joint (B, CvL) templates, and runs the closure test.  The resulting
# templates.root is what cards/belphes_card_CMS_allflavors.tcl points at.
#
# Requires ROOT on PATH (root-config).  Under CMSSW:  cmsenv
#
# Usage:
#   ./build_templates.sh --fetch 4 --smallest      # quick end-to-end test
#   ./build_templates.sh --fetch 68                # full 17.6 GB record
#   ./build_templates.sh                           # reuse whatever is in data/
#   ./build_templates.sh --config coarse.conf -o coarse.root
#
set -euo pipefail
cd "$(dirname "$0")"

CONFIG="binning.conf"
OUTPUT="templates.root"
FETCH=""
SMALLEST=""
NOVALIDATE=0

usage() { sed -n '3,16p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--config)  CONFIG="$2"; shift 2 ;;
    -o|--output)  OUTPUT="$2"; shift 2 ;;
    --fetch)      FETCH="$2";  shift 2 ;;
    --smallest)   SMALLEST="--smallest"; shift ;;
    --no-validate) NOVALIDATE=1; shift ;;
    -h|--help)    usage 0 ;;
    *) echo "unknown option: $1" >&2; usage 1 ;;
  esac
done

command -v root-config >/dev/null || {
  echo "ERROR: root-config not found. Set up ROOT first (under CMSSW: cmsenv)." >&2
  exit 1
}

# ---- compile ---------------------------------------------------------------
# Rebuild only when the source is newer than the binary.
for tool in build_templates validate_templates; do
  if [[ ! -x "$tool" || "$tool.cpp" -nt "$tool" ]]; then
    echo "==> compiling $tool"
    g++ -O2 -Wall "$tool.cpp" $(root-config --cflags --libs) -o "$tool"
  fi
done

# ---- step 1: fetch ---------------------------------------------------------
if [[ -n "$FETCH" ]]; then
  ./fetch_opendata.sh --limit "$FETCH" $SMALLEST
fi

shopt -s nullglob
INPUTS=(data/*.root)
shopt -u nullglob
if [[ ${#INPUTS[@]} -eq 0 ]]; then
  echo "ERROR: no .root files in data/. Run with --fetch N first." >&2
  exit 1
fi

# ---- step 2: build ---------------------------------------------------------
echo
echo "==> building templates from ${#INPUTS[@]} file(s)"
./build_templates -c "$CONFIG" -o "$OUTPUT" "${INPUTS[@]}"

# ---- closure test ----------------------------------------------------------
if [[ "$NOVALIDATE" -eq 0 ]]; then
  for flavor in 5 4 0; do
    echo
    echo "==> closure test, hadronFlavour $flavor"
    ./validate_templates -t "$OUTPUT" -f "$flavor" "${INPUTS[@]}"
  done
fi

echo
echo "==> $OUTPUT ready. Point TemplateFile at it in your Delphes card."
