#!/usr/bin/env bash
# Fetch raw third-party archives required for full local pipeline reproduction:
#   1. Open Khipu Repository (OKR v2.1.0, Zenodo DOI 10.5281/zenodo.18025748) -> data/okr/
#   2. Khipu Field Guide Excel workbooks (khipufieldguide/data-science-khipu-article) -> data/kfg/
#   3. Optional (--with-textos): raw Textos Andinos recitation pages -> data/textos_andinos_raw/
#      (Not redistributed in Git; subject to upstream copyright/redistribution terms.)

set -euo pipefail

mkdir -p data/okr data/kfg

echo "[1/2] Downloading Open Khipu Repository (OKR v2.1.0) from Zenodo..."
curl -fL "https://zenodo.org/api/records/18025748/files-archive" -o data/okr_v210.zip
unzip -qo data/okr_v210.zip -d data/okr/
# Handle nested release archive if present
for z in data/okr/*.zip; do
  [ -f "$z" ] && unzip -qo "$z" -d data/okr/
done
rm -f data/okr_v210.zip

echo "[2/2] Downloading Khipu Field Guide (KFG) Excel workbooks from GitHub..."
curl -fL "https://github.com/khipufieldguide/data-science-khipu-article/archive/refs/heads/main.zip" -o data/kfg_main.zip
unzip -qo data/kfg_main.zip -d data/kfg_tmp/
if [ -d "data/kfg_tmp/data-science-khipu-article-main/data" ]; then
  cp -r data/kfg_tmp/data-science-khipu-article-main/data data/kfg/
fi
rm -rf data/kfg_tmp data/kfg_main.zip

if [ "${1:-}" = "--with-textos" ]; then
  echo "[Optional] Fetching Textos Andinos recitation pages into data/textos_andinos_raw/..."
  mkdir -p data/textos_andinos_raw
  for id in m_04 m_08 m_09 m_10 m_14 m_15 m_16 m_23 m_24 m_25 m_26 m_48 m_60 m_63; do
    curl -fsSL "https://www.khipufieldguide.com/textosbook/recitations/${id}.html" -o "data/textos_andinos_raw/${id}.html" || true
  done
fi

echo "Done. Verify OKR SQLite DB under data/okr/ and KFG workbooks under data/kfg/data/KFG_DB_Excel/."
