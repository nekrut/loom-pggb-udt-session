# Loom session: PGGB pangenome analysis as a Galaxy user-defined tool

A Loom/Orbit working directory that drove a `pggb` (Pangenome Graph
Builder) workflow against `test.galaxyproject.org` as a **user-defined
tool (UDT)**, with jobs dispatched through the Galaxy relay
(`relay.usegalaxy.org`) to a **local `pulsar` endpoint** running on the
analyst's machine. Tool processes executed inside Docker on that
machine; the `pggb`, `odgi`, and synteny-plotting steps therefore ran
locally while staying inside Galaxy's history and provenance model.

## What this session produced

We registered nine user-defined tools in the Galaxy account associated
with this session — a single-column cutter, a `pggb` aligner, an
`odgi viz` per-chromosome wrapper iterated to v0.4, and a
GFA → synteny-strip plotter. Each iteration created a new UDT version
(Galaxy has no in-place update API) under five distinct `tool_id`s:
`cut_column`, `pggb_align`, `odgi_viz`, `odgi_viz_chr`,
`synteny_from_gfa`. The accumulated UDT versions and their UUIDs are
recoverable from the Loom session log
(`activity.jsonl`, `session.jsonl`).

The first UDT-creation call returned `HTTP 403 "User is not allowed to
run unprivileged"` — the canonical signature of TPV rejecting
`tool_type_user_defined`. Subsequent calls succeeded once the
admin-side destination accepted the user-defined tool type. This is
captured upstream in
[galaxyproject/loom#98](https://github.com/galaxyproject/loom/issues/98)
as a setup gotcha that the seamless first-run flow must surface.

## Method

We simulated five iterations of small synthetic genomes diverged from
a reference at controlled rates (1, 2, 5, 10, 15%), built pangenome
graphs with `pggb`, and compared the resulting graph and synteny
structure against the ground-truth divergence parameters.

| Iteration | Directory | Notes |
|-----------|-----------|-------|
| 1 | `genomes/` | Initial `pangenome.fa`, no graph yet. |
| 2 | `genomes2/` | `pangenome2.fa`. |
| 3 | `genomes3/` | First with truth files: `truth.tsv`, `graph_identity.tsv`, `synteny_pairwise_truth.png`, `synteny_truth.png`, `pangenome3.gfa`, `identity_strip.png`. |
| 4 | `genomes4/` | Larger / parameter sweep. |
| 5 | `genomes5/` | Final iteration. |

Per-divergence read alignments and variant calls land in
`igv/{div01,div02,div05,div10,div15}.{bam,vcf.gz}` (with `.bai` /
`.tbi` indices) for inspection in IGV.

## Layout

```
├── notebook.md                — Loom durable project log (empty in this snapshot).
├── echydna.py                 — Top-level driver / orchestration helper.
├── simulate_genomes*.py       — Synthetic-genome generators (base, indels, rearrangements, 1Mb variants).
├── graph_identity.py          — Identity statistics from a pangenome graph.
├── plot_strips.py             — Identity / divergence strip plot helper.
├── synteny_pairwise*.py       — Pairwise synteny plot variants (incl. SynVisio-style).
├── synteny_plot.py            — Synteny visualization driver.
├── synteny_from_gfa.py        — Strip plot from a GFA pangenome graph (mirrors the UDT of the same name).
├── regen_with_truth.py        — Regenerate analysis bundling the ground-truth comparison.
├── end2end_*.gfa, end2end_synteny_*.png
│                              — End-to-end pangenome + synteny outputs across runs.
├── pangenome*_chr*.png, galaxy_synteny_pangenome*.png
│                              — Chromosome-level visualisations from intermediate iterations.
├── genomes / genomes2…5 /     — Per-iteration inputs and outputs (table above).
├── igv/                       — BAM + VCF tracks per divergence rate.
├── activity.jsonl             — Loom activity log for this session (601 events).
└── session.jsonl              — Symlink to the pi-agent session log under `~/.pi/agent/sessions/…`.
                                 (Broken on a clone; preserved for reference.)
```

## Reproducing the runtime

The execution path was:

```
Loom / Orbit  →  test.galaxyproject.org  →  relay.usegalaxy.org
                                              ↓
                                   Local pulsar (manager 'anton')
                                              ↓
                                       Docker on this machine
```

Setup notes captured in the upstream issues:

- [`galaxyproject/loom#97`](https://github.com/galaxyproject/loom/issues/97)
  — broad capability + upstream `pulsar` blockers
  ([`galaxyproject/pulsar#452`](https://github.com/galaxyproject/pulsar/issues/452)).
- [`galaxyproject/loom#98`](https://github.com/galaxyproject/loom/issues/98)
  — client-side first-run setup for the local-Pulsar mode (auto-config
  of `app.yml` from Galaxy preferences, supervised local
  `pulsar-main`, Docker-group handling).
- [`galaxyproject/galaxy-skills#17`](https://github.com/galaxyproject/galaxy-skills/pull/17)
  — `user-defined-tools` skill nudging a minimal `help` section on
  every UDT (motivated directly by this session, where none of the
  registered UDTs included one).

A non-trivial finding: `command.sh` injected by Galaxy hard-codes
`GALAXY_SLOTS=1` when no recognised scheduler is detected (no Slurm,
SGE, PBS, LSF, k8s, HTCondor). Because Pulsar runs `command.sh` as a
plain bash process, `--cpus ${GALAXY_SLOTS:-1}` always evaluates to
`--cpus 1`, capping every container at one core regardless of host
capacity. Lifting this requires both a destination-side `cores` bump
*and* an injection of `GALAXY_SLOTS` into the Pulsar job environment.

## Scripts of note

| Script | Role |
|--------|------|
| `simulate_genomes.py`, `simulate_genomes_indels.py`, `simulate_genomes_indels_1mb.py`, `simulate_genomes_rearr.py` | Generate synthetic genome sets at controlled divergence rates / variant types. |
| `graph_identity.py` | Compute graph-level identity statistics from a `pggb` GFA. |
| `synteny_from_gfa.py` | GFA → synteny-strip plot; mirrors the `synteny_from_gfa` UDT registered in Galaxy. |
| `synteny_pairwise*.py`, `synteny_plot.py` | Pairwise synteny visualisations, including a SynVisio-style variant. |
| `regen_with_truth.py` | Regenerate analysis with the ground-truth divergence comparison overlaid. |
| `echydna.py` | Driver / orchestration helper threading the simulation → `pggb` → analysis pipeline. |

## Status

Snapshot of an exploratory session. `notebook.md` is empty in this
revision; `activity.jsonl` carries the timestamped tool-call record
that drove the work.
