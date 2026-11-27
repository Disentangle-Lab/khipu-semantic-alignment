"""
build_colonial_schemas.py — Reproducibly writes the primary-source verified colonial
khipu transcription schemas (JSON) into data/colonial_schemas/.

Every cell is transcribed verbatim from Textos Andinos (Pärssinen & Kiviharju;
Manuel Medrano's 72-document corpus hosted at https://www.khipufieldguide.com/textosbook/).
No columns are cherry-picked or imputed; editorial discrepancies are explicitly documented.
"""

from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "colonial_schemas"
OUT.mkdir(parents=True, exist_ok=True)

SCHEMAS = [
    {
        "schema_id": "chucuito_1567_moiety_k5",
        "title": "Chucuito 4-Khipu Tributary & Population Account (Martín Cari, 1532/1567)",
        "source_citation": "Textos Andinos m_23, m_24, m_25, m_26 (Visita de Chucuito, Garci Diez de San Miguel 1567; KFG textosbook)",
        "region": "Chucuito (Lake Titicaca / Collao)",
        "year": 1567,
        "domain": "census_tribute",
        "role": "step3_primary",
        "editorial_notes": (
            "12 rows = 7 cabeceras split by Anansaya/Lurinsaya moiety (Ilave and Yunguyo are reported combined "
            "across both moieties in m_23, m_24, m_26; m_25 splits them 792+792=1584 and 561+594=1155). "
            "Juli Anansaya Aymaras includes 153 Chinchasuyo mitimaes (1438+153=1591) and Minas includes 3 mitimaes (45+3=48). "
            "Pomata Anansaya Aymaras includes 20 Canas mitimaes (1663+20=1683). Sama (yungas) omitted from 12-row table "
            "because it lacks m_24/m_26 entries."
        ),
        "column_names": ["aymaras", "uros", "ropa_piezas", "tributo_pesos", "indios_minas"],
        "row_names": [
            "Chucuito_Anansaya", "Chucuito_Lurinsaya",
            "Acora_Anansaya", "Acora_Lurinsaya",
            "Ilave_Ambas",
            "Juli_Anansaya", "Juli_Lurinsaya",
            "Pomata_Anansaya", "Pomata_Lurinsaya",
            "Yunguyo_Ambas",
            "Zepita_Anansaya", "Zepita_Lurinsaya",
        ],
        "value_matrix": [
            [1233, 500, 83, 1418, 41],
            [1384, 347, 83, 1418, 41],
            [1221, 440, 77, 1320, 39],
            [1207, 378, 76, 1320, 38],
            [1470, 1070, 93, 1584, 46],
            [1591, 158, 92, 1518, 48],
            [1804, 256, 114, 1814, 57],
            [1683, 110, 106, 1749, 53],
            [1341, 183, 85, 1452, 42],
            [1039, 381, 66, 1155, 33],
            [1112, 186, 71, 1188, 35],
            [866, 120, 54, 924, 27],
        ],
        "structural_flags": {
            "has_moiety_split": True,
            "has_summary_row": False,
            "has_round_100_totals": False,
        },
    },
    {
        "schema_id": "huaraz_1558_demo_k6",
        "title": "Visita de Huaraz Demographic Census by Principal (1558, 6 person categories)",
        "source_citation": "Textos Andinos m_63 (Visita de los repartimientos de Guaraz y Llaguaraz, Lic. Diego Álvarez 1558; KFG textosbook)",
        "region": "Huaraz, Ancash (Callejón de Huaylas / Santa Basin)",
        "year": 1558,
        "domain": "census",
        "role": "step3_primary",
        "editorial_notes": (
            "Verbatim counts for all 7 principales recorded in m_63, taking the 6 strictly demographic person categories "
            "in verbatim recitation order (excluding casas). Martín Yarollivia reports no mozos (0). "
            "Document ends with an explicit repartimiento-wide summary block."
        ),
        "column_names": ["casados", "mozos", "mozas", "muchachos_6_13", "ninos_0_5", "viejos_viudos_50plus"],
        "row_names": [
            "Don_Diego_Chinchecaque",
            "Santiago_Xulcallivia",
            "Martin_Yarollivia",
            "Geronimo_Rimaycochache",
            "Don_Gonzalo_Xulcacollas",
            "Domingo_Vilcacaque",
            "Don_Fernando_Paverr",
        ],
        "value_matrix": [
            [27, 2, 12, 28, 24, 14],
            [41, 4, 11, 30, 34, 22],
            [34, 0, 7, 21, 18, 22],
            [39, 3, 5, 35, 33, 23],
            [32, 2, 8, 26, 24, 24],
            [27, 3, 11, 21, 17, 21],
            [43, 5, 8, 26, 24, 11],
        ],
        "structural_flags": {
            "has_moiety_split": False,
            "has_summary_row": True,
            "has_round_100_totals": False,
        },
    },
    {
        "schema_id": "huaraz_1558_full_k7",
        "title": "Visita de Huaraz Full Census by Principal (1558, casas + 6 person categories)",
        "source_citation": "Textos Andinos m_63 (Visita de los repartimientos de Guaraz y Llaguaraz, Lic. Diego Álvarez 1558; KFG textosbook)",
        "region": "Huaraz, Ancash (Callejón de Huaylas / Santa Basin)",
        "year": 1558,
        "domain": "census",
        "role": "step3_primary",
        "editorial_notes": "Identical to huaraz_1558_demo_k6 with the leading casas column preserved verbatim.",
        "column_names": ["casas", "casados", "mozos", "mozas", "muchachos_6_13", "ninos_0_5", "viejos_viudos_50plus"],
        "row_names": [
            "Don_Diego_Chinchecaque",
            "Santiago_Xulcallivia",
            "Martin_Yarollivia",
            "Geronimo_Rimaycochache",
            "Don_Gonzalo_Xulcacollas",
            "Domingo_Vilcacaque",
            "Don_Fernando_Paverr",
        ],
        "value_matrix": [
            [44, 27, 2, 12, 28, 24, 14],
            [58, 41, 4, 11, 30, 34, 22],
            [44, 34, 0, 7, 21, 18, 22],
            [65, 39, 3, 5, 35, 33, 23],
            [56, 32, 2, 8, 26, 24, 24],
            [47, 27, 3, 11, 21, 17, 21],
            [71, 43, 5, 8, 26, 24, 11],
        ],
        "structural_flags": {
            "has_moiety_split": False,
            "has_summary_row": True,
            "has_round_100_totals": False,
        },
    },
    {
        "schema_id": "huaraz_1532_full_k7",
        "title": "Visita de Huaraz Pre-Conquest / Early Census by Principal (m_60, 7 complete principales)",
        "source_citation": "Textos Andinos m_60 (Visita de los repartimientos de Guaraz y Llaguaraz, Lic. Diego Álvarez; KFG textosbook)",
        "region": "Huaraz, Ancash (Callejón de Huaylas / Santa Basin)",
        "year": 1532,
        "domain": "census",
        "role": "step3_primary",
        "editorial_notes": (
            "Verbatim 7-column counts for the 7 principales with complete declarations in m_60. "
            "Principal #1 (Don Hernando Guamacochache) lacks the initial casas/casados header line in the surviving text "
            "and is excluded rather than imputed; Pedro Yana (mitimaes Yauyos, 8 casas, 5 casados, 3 mozos, 0 mozas, "
            "4 muchachos, 14 ninos, 5 viudos) is omitted from the 7-row core because mozas is unmentioned."
        ),
        "column_names": ["casas", "casados", "mozos", "mozas", "muchachos_6_13", "ninos_0_5", "viejos_viudos_50plus"],
        "row_names": [
            "Don_Martin_Guache",
            "Santiago_Pomaguaman",
            "Don_Antonio_Paucarllangue",
            "Don_Gonzalo_Yaroyana",
            "Don_Diego_Rahe",
            "Don_Pedro_Ricari",
            "Don_Felipe_Caxayaure",
        ],
        "value_matrix": [
            [45, 32, 5, 12, 27, 56, 27],
            [88, 59, 6, 22, 35, 86, 33],
            [37, 42, 11, 8, 22, 65, 19],
            [48, 44, 7, 13, 23, 49, 14],
            [31, 25, 2, 6, 14, 25, 5],
            [57, 36, 3, 16, 31, 80, 31],
            [69, 46, 5, 9, 11, 55, 22],
        ],
        "structural_flags": {
            "has_moiety_split": False,
            "has_summary_row": False,
            "has_round_100_totals": False,
        },
    },
    {
        "schema_id": "sacaca_1551_k22",
        "title": "Sacaca 4-Year Encomienda Tribute Account (m_14 / m_15 / m_16 consensus, 1551)",
        "source_citation": "Textos Andinos m_14 (Hernando Achacata), m_15 (Luis Anba), m_16 (Pedro Horuro), Sacaca 1551; KFG textosbook",
        "region": "Sacaca (Charcas / Potosí)",
        "year": 1551,
        "domain": "storehouse_tribute",
        "role": "step3_primary",
        "editorial_notes": (
            "Verbatim 4-year x 22-item canonical sequence recited identically by three independent khipukamayuqs "
            "(m_14, m_15, m_16). Uses majority/consensus across the three witnesses: yr-1 silver=21200 (m_14 & m_16; m_15=28200), "
            "yr-1 ovejas=500 (m_15 & m_16; m_14=300), yr-1 gallinas=103 (m_15 & m_16; m_14=130), yr-2 maiz=550 (m_15 & m_16; m_14=250), "
            "yr-2 sebo=16 (explicit in m_16; omitted in m_14/m_15), yr-2 huevos=1304 (m_15 & m_16; m_14=1004), "
            "yr-3 quinua=1.5 fanegas / 3 cargas (m_15=3 cargas; m_14/m_16=1.5 fanegas -> 3), yr-4 cera=4 (m_15; m_14=2, m_16=12), "
            "yr-4 sal=6 (m_14 & m_16; m_15=5). All 22 columns appear in fixed canonical order every year."
        ),
        "column_names": [
            "pesos_plata", "indios_minas", "maiz_fanegas", "chuno_fanegas",
            "ovejas_carneros", "carneros_comer", "corderos", "sebo_arrobas",
            "ropa_abasca", "frazadas", "lana_arrobas", "mandiles",
            "jaquimas_cinchas", "pescuezos", "miel_arrobas", "cera_arrobas",
            "gallinas", "patos", "huevos", "quinua_cargas", "papas_cargas", "sal_cargas",
        ],
        "row_names": ["Year_1", "Year_2", "Year_3", "Year_4"],
        "value_matrix": [
            [21200, 100, 750, 250, 500, 52, 53, 20, 120, 60, 10, 50, 80, 40, 20, 23, 103, 208, 1205, 15, 10, 30],
            [16960, 80, 550, 200, 350, 52, 53, 16, 90, 35, 7, 30, 45, 21, 14, 8, 60, 130, 1304, 10, 13, 23],
            [14840, 70, 450, 125, 150, 42, 31, 7, 35, 21, 3, 11, 13, 6, 7, 4, 30, 50, 850, 3, 4, 11],
            [12432, 60, 300, 100, 100, 21, 19, 3, 20, 11, 1, 6, 11, 4, 4, 4, 15, 34, 540, 1, 5, 6],
        ],
        "structural_flags": {
            "has_moiety_split": False,
            "has_summary_row": False,
            "has_round_100_totals": True,
        },
    },
]


def main():
    for s in SCHEMAS:
        path = OUT / f"{s['schema_id']}.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(s, f, indent=2, ensure_ascii=False)
        mat = s["value_matrix"]
        print(f"Wrote {path.name:32s} ({len(mat):2d} rows x {len(mat[0]):2d} cols)")


if __name__ == "__main__":
    main()
