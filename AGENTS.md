# Developer notes

Read `docs/adr/` before any architectural change. This file is a condensed index; the ADRs are the source of truth.

## Architecture (ADR-0001)

Three layers, strict inward dependency:
- `src/application/` — use cases, orchestration, CLI
- `src/domain/` — pure business logic, no external dependencies
- `src/shared/` — cross-cutting utilities, exceptions, logger

Business logic always lives in Domain. Application only orchestrates.

## Agentic cell pattern (ADR-0002)

Each cell is an autonomous agent with two accumulators:
- **D1** (DNA damage): threshold 2.0 → PRIMER stage
- **D2** (immunosuppression): threshold 5.0 → TUMORAL stage

Five phases in `live()`: baseline check → intrinsic apoptosis (TP53-based) → extrinsic apoptosis (D2-based) → neoplasm development → division/remodeling.

Six `CellLifeStage` values: `DEAD, BASELINE, UNSTABLE, UNPROTECTED, PRIMER, TUMORAL`.

**Critical rule:** BRCA1's native state is always `+/-` (heterozygous), never `+/+`. It is a tumor suppressor; the native carrier state has one copy already lost.

## Mutations (ADR-0003)

- TP53: `+/+` → `+/-` → `-/-` (bidirectional progression toward loss of function)
- BRCA1: `+/-` → `-/-` (unidirectional)
- Mutation probability: `threshold + instability_k * genomic_instability`
- D1/D2 are inherited by daughter cells on division

## Apoptosis (ADR-0004)

- **Intrinsic** (TP53-based): triggered in Phase 2, probability scales with TP53 loss
- **Extrinsic** (D2-based): triggered in Phase 3; below the D2 threshold the cell is visible to immune clearance, above it the cell is effectively camouflaged
- **Lethality rule:** BRCA1 -/- is lethal only if TP53 is +/+ or +/-. It is NOT lethal if TP53 is already -/- (TP53 status controls lethality, not BRCA1 alone)

## Big Bang mode (ADR-0005)

Neoplastic cells divide at `neoplastic_division_rate` instead of the baseline `division_rate` once transformed, modeling explosive clonal expansion. Governed by `enable_big_bang_mode`. See `docs/CALIBRATION_AND_PHENOTYPES.md` for why this mechanism is structurally necessary to reproduce clinical incidence curves.

## Tissue manager (ADR-0006)

`Tissue` owns cells via `unique_ptr`, runs the lifecycle loop, and catches cell-death exceptions and division signals to add/remove cells safely.

## Rules

| Do | Don't |
|---|---|
| Business logic in Domain | Business logic in Application |
| `unique_ptr` for ownership | Raw pointers (except in getters) |
| BRCA1 native state = `+/-` | BRCA1 native state = `+/+` |
| Exceptions for cell death | Silent state mutation on death |
| `ICell` interface, not concrete types | Casting to a concrete cell class outside tests |

## Where to look

| Question | Look at |
|---|---|
| Why did a cell die with no obvious cause? | Phase 1/2 apoptosis, TP53 status, `alive()` |
| Why isn't a cell dividing? | D1 < 2.0, or TP53 not yet -/- |
| Why is a cell dividing very fast? | Big Bang mode + `is_neoplastic_` |
| Why is D1 growing fast? | BRCA1 -/-, high `instability_k` |
| Why won't a cell reach PRIMER? | D1 still ≤ 2.0 even with TP53 -/- |
| Why does a neoplastic cell still die? | D2 < 5.0, still subject to extrinsic apoptosis |

## Known documentation gap

`docs/README.md` describes an earlier model revision (e.g. instant BRCA1 -/- lethality in a fixed phase) that conflicts with the calibrated model described in the top-level `README.md` and `docs/CALIBRATION_AND_PHENOTYPES.md`. Needs reconciling — see the note in the top-level README.
