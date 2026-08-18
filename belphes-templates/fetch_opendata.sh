#!/usr/bin/env bash
#
# fetch_opendata.sh -- step 1 of the belphes template pipeline.
#
# Queries the CERN Open Data API for a record, resolves the NanoAOD file list,
# and downloads files over XRootD.  The full TTJets record (67727) is 17.6 GB
# across 68 files, so --limit / --smallest exist to pull a testable subset.
#
# Usage:
#   ./fetch_opendata.sh --list
#   ./fetch_opendata.sh --limit 1 --smallest        # ~5.8 MB, good for testing
#   ./fetch_opendata.sh --limit 10 --outdir data
#   ./fetch_opendata.sh --all                       # the whole 17.6 GB
#
set -euo pipefail

RECORD=67727
OUTDIR="$(dirname "$0")/data"
LIMIT=""
SMALLEST=0
LIST_ONLY=0
FETCH_ALL=0

usage() { sed -n '3,15p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --record)   RECORD="$2"; shift 2 ;;
    --outdir)   OUTDIR="$2"; shift 2 ;;
    --limit)    LIMIT="$2";  shift 2 ;;
    --smallest) SMALLEST=1;  shift ;;
    --list)     LIST_ONLY=1; shift ;;
    --all)      FETCH_ALL=1; shift ;;
    -h|--help)  usage 0 ;;
    *) echo "unknown option: $1" >&2; usage 1 ;;
  esac
done

if [[ -z "$LIMIT" && "$FETCH_ALL" -eq 0 && "$LIST_ONLY" -eq 0 ]]; then
  echo "refusing to download without --limit N or --all (record is O(10) GB)" >&2
  exit 1
fi

command -v xrdcp   >/dev/null || { echo "xrdcp not found (need XRootD client)" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found" >&2; exit 1; }

mkdir -p "$OUTDIR"
META="$OUTDIR/.record_${RECORD}.json"

if [[ ! -s "$META" ]]; then
  echo "==> fetching metadata for record $RECORD"
  curl -sSf --max-time 120 "https://opendata.cern.ch/api/records/${RECORD}" -o "$META"
fi

# size<TAB>filename<TAB>uri, one line per file
INDEX="$OUTDIR/.files_${RECORD}.tsv"
python3 - "$META" > "$INDEX" <<'PYEOF'
import json, sys
meta = json.load(open(sys.argv[1]))["metadata"]
files = [f for idx in meta.get("_file_indices", []) for f in idx.get("files", [])]
files += [f for f in meta.get("files", []) if f.get("uri", "").endswith(".root")]
seen = set()
for f in files:
    uri = f.get("uri", "")
    if not uri.endswith(".root") or uri in seen:
        continue
    seen.add(uri)
    print("%d\t%s\t%s" % (f.get("size", 0), f.get("filename", uri.rsplit("/", 1)[-1]), uri))
PYEOF

NFILES=$(wc -l < "$INDEX")
TOTAL=$(awk -F'\t' '{s+=$1} END {printf "%.1f", s/1e9}' "$INDEX")
echo "==> record $RECORD: $NFILES ROOT files, ${TOTAL} GB total"
[[ "$NFILES" -gt 0 ]] || { echo "no ROOT files found in record metadata" >&2; exit 1; }

# Selection order: by size ascending if --smallest, else as published.
SELECTED="$OUTDIR/.selected_${RECORD}.tsv"
if [[ "$SMALLEST" -eq 1 ]]; then sort -n -k1,1 "$INDEX" > "$SELECTED"; else cp "$INDEX" "$SELECTED"; fi
if [[ -n "$LIMIT" ]]; then head -n "$LIMIT" "$SELECTED" > "$SELECTED.tmp" && mv "$SELECTED.tmp" "$SELECTED"; fi

if [[ "$LIST_ONLY" -eq 1 ]]; then
  awk -F'\t' '{printf "  %8.1f MB  %s\n", $1/1e6, $2}' "$SELECTED"
  exit 0
fi

NSEL=$(wc -l < "$SELECTED")
SELGB=$(awk -F'\t' '{s+=$1} END {printf "%.2f", s/1e9}' "$SELECTED")
echo "==> downloading $NSEL file(s), ${SELGB} GB -> $OUTDIR"

ok=0; skip=0; fail=0
while IFS=$'\t' read -r size name uri; do
  dest="$OUTDIR/$name"
  # Resume-safe: treat a file of the published size as already done.
  if [[ -f "$dest" && "$(stat -c%s "$dest")" == "$size" ]]; then
    echo "  [skip] $name"; skip=$((skip+1)); continue
  fi
  echo "  [get ] $name ($(awk -v s="$size" 'BEGIN{printf "%.1f", s/1e6}') MB)"
  if xrdcp -f -s "$uri" "$dest"; then ok=$((ok+1)); else
    echo "  [FAIL] $name" >&2; rm -f "$dest"; fail=$((fail+1))
  fi
done < "$SELECTED"

echo "==> done: $ok downloaded, $skip already present, $fail failed"
[[ "$fail" -eq 0 ]]
