#!/usr/bin/env bash
# Fetch raw third-party archives required for full local pipeline reproduction
# with pinned versions and SHA-256 integrity verification:
#   1. Open Khipu Repository (OKR v2.1.0, Zenodo DOI 10.5281/zenodo.18025748) -> data/okr/
#   2. Khipu Field Guide Excel workbooks (pinned commit 167b8df031dc9b8f1819c7a1148e35836d02016b) -> data/kfg/
#   3. Optional (--with-textos): raw Textos Andinos recitation pages -> data/textos_andinos_raw/
#      (Not redistributed in Git; subject to upstream copyright/redistribution terms.)

set -euo pipefail

OKR_URL="https://zenodo.org/api/records/18025748/files/open-khipu-repository-2.1.0.zip/content"
OKR_SHA256="66cad1ec007ca00649326255b3ccb74a6459bf36e50bff9b08b74fc1ea09402e"

KFG_COMMIT="167b8df031dc9b8f1819c7a1148e35836d02016b"
KFG_URL="https://github.com/khipufieldguide/data-science-khipu-article/archive/${KFG_COMMIT}.zip"
KFG_SHA256="fa58f03c7be5b8eb1691f4506e9376e7df0b40e1f1a812a16cdd98e32b4d67e9"

verify_sha256() {
  local file="$1"
  local expected="$2"
  local actual
  if command -v sha256sum >/dev/null 2>&1; then
    actual="$(sha256sum "$file" | awk '{print $1}')"
  elif command -v shasum >/dev/null 2>&1; then
    actual="$(shasum -a 256 "$file" | awk '{print $1}')"
  else
    echo "ERROR: Neither sha256sum nor shasum found to verify $file" >&2
    exit 1
  fi
  if [ "$actual" != "$expected" ]; then
    echo "ERROR: SHA-256 mismatch for $file" >&2
    echo "  expected: $expected" >&2
    echo "  actual:   $actual" >&2
    exit 1
  fi
  echo "  SHA-256 verified OK ($actual)"
}

mkdir -p data/okr data/kfg

echo "[1/2] Downloading Open Khipu Repository (OKR v2.1.0) from Zenodo..."
curl -fL "$OKR_URL" -o data/okr_v210.zip
verify_sha256 data/okr_v210.zip "$OKR_SHA256"
unzip -qo data/okr_v210.zip -d data/okr/
# Handle nested release archive if present
for z in data/okr/*.zip; do
  [ -f "$z" ] && unzip -qo "$z" -d data/okr/
done
rm -f data/okr_v210.zip

echo "[2/2] Downloading Khipu Field Guide (KFG) at pinned commit ${KFG_COMMIT}..."
curl -fL "$KFG_URL" -o data/kfg_pinned.zip
verify_sha256 data/kfg_pinned.zip "$KFG_SHA256"
unzip -qo data/kfg_pinned.zip -d data/kfg_tmp/
if [ -d "data/kfg_tmp/data-science-khipu-article-${KFG_COMMIT}/data" ]; then
  cp -r "data/kfg_tmp/data-science-khipu-article-${KFG_COMMIT}/data" data/kfg/
fi
rm -rf data/kfg_tmp data/kfg_pinned.zip

if [ "${1:-}" = "--with-textos" ]; then
  echo "[Optional] Fetching Textos Andinos recitation pages into data/textos_andinos_raw/..."
  mkdir -p data/textos_andinos_raw
  for id in m_04 m_08 m_09 m_10 m_14 m_15 m_16 m_23 m_24 m_25 m_26 m_48 m_60 m_63; do
    curl -fsSL "https://www.khipufieldguide.com/textosbook/recitations/${id}.html" -o "data/textos_andinos_raw/${id}.html" || true
  done
fi

echo "Done. Verify OKR SQLite DB under data/okr/ and KFG workbooks under data/kfg/data/KFG_DB_Excel/."
