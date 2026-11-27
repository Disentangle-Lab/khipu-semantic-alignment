# Reproducible pipeline Makefile (fixed deterministic seeds inside each module:
# Step 1/Hardening SEED=2026, Step 3 SEED=1670, Phase 2A/2C SEED=42, P0 Audits SEED=20261127).

.PHONY: all prereg phase1 phase2 figures test verify-hash

verify-hash:
	sha256sum -c PREREGISTRATION.sha256

test:
	uv run pytest -q

phase1:
	uv run python -m src.step1_archive_invariance
	uv run python -m src.step2_taxonomy_outliers
	uv run python -m src.step1_2_hardening
	uv run python -m src.gate0_santa_valley
	uv run python -m src.build_colonial_schemas
	uv run python -m src.step3_colonial_alignment

phase2:
	uv run python -m src.phase2a_subsidiary_grammar
	uv run python -m src.phase2c_pair_discovery
	uv run python -m src.p0_audits_and_power

figures:
	uv run --group figures python -m src.make_figures

all: verify-hash phase1 figures phase2 test
