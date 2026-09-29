# Katz-BP decision-impact pilot

This pilot answers one question: **does the bounded Katz-centrality term change which patient
class gets served, compared with plain backpressure (BP)?** It does **not** test whether
congestion improves.

Registration and triage are served FCFS (`model.fcfs_nodes`). The policy decides only at
fast_track, consultation, lab, imaging and treatment, and only those nodes produce decision
epochs.

The ED is simulated under PlainBP. At every decision epoch, the KatzBP choice for each
(κ, θ) is computed on the same live state without acting on it (shadow evaluation). The
pilot then measures how often the two choices differ.

## How to run

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m pytest            # all tests must pass first
python -m src.pilot --config configs/pilot.yaml   # writes to results/
```

`--out DIR` sends the outputs somewhere other than `experiment.output_dir`. A full run
(5 replications × 8 simulated days × 43 shadow policies: 36 KatzBP (3 κ × 12 θ) + 7 OneHopBP)
takes about 5 s on a laptop. The scenario comparison adds another 4 scenarios × 5 replications, so the whole run
takes about 10 s.

## Layout

| Path | Contents |
|---|---|
| `configs/pilot.yaml` | Every parameter. Model numbers are marked `PROVISIONAL — replace after calibration`. |
| `src/model.py` | Nodes, classes, per-stage buffers, class mix, and analytic λ calibration. |
| `src/network.py` | Transition-frequency adjacency A, ρ(A), the Katz solve with the bounded ĉ, and the one-hop score ĥ. |
| `src/policies.py` | The policy interface `select(node, state, rng) -> class_id`, plus `FCFS` (used at the FCFS nodes), `PlainBP`, `KatzBP(theta, kappa)` and `OneHopBP(theta_h, kappa)`. |
| `src/sim.py` | The SimPy event-driven simulation, RNG streams, the decision-epoch hook, blocking-after-service with finite buffers, and the service-rate shock. |
| `src/shadow.py` | The shadow evaluator that writes the per-epoch log. |
| `src/metrics.py` | Per-replication metrics, mean ± 95% CI, and the automatic verdict. |
| `src/report.py`, `src/plots.py` | The markdown report, the console summary, and the figures. |
| `src/scenarios.py` | The scenario comparison (A–D and C-window): config merging, runs, metrics and per-scenario verdicts. |
| `src/pilot.py` | The CLI and the orchestration. |
| `tests/` | The pytest suite (spec §8). |

Adding APQ later means subclassing `Policy` (or `WeightedPolicy` for an
argmax rule) in `src/policies.py`. You then pass the new policy to `EDSimulation` as the
acting policy, or to `ShadowLogger` as a shadow.

## Outputs (`results/`)

| File | Contents |
|---|---|
| `pilot_report.md` | Verdict at the top (including the per-scenario verdicts), then the candidate-set table, set-up, all metric tables, the figures, "Katz vs one-hop" and "Scenario comparison". |
| `pilot_summary.csv` | Long format with one row per (metric, κ, θ, node): `mean, ci_halfwidth, ci_low, ci_high, n_reps`. `node = ALL` means overall; an empty κ or θ means the metric does not depend on it. |
| `epochs/epochs_rep{r}.csv.gz` | The post-warm-up per-epoch log (spec §5). |
| `scenario_summary.csv` | The scenario comparison in long format: `scenario` (A, B, C, C-window, D), `metric`, `node`, `mean`, CI and `n_reps`. |
| `scenarios/{A,B,C,D}/epochs_rep{r}.csv.gz` | The per-epoch logs for each scenario (C-window is a time filter of C). |
| `classes.csv` | Maps `class_id` (buffer id) to level, pathway, stage, node, next node and mean service time. |
| `resolved_config.yaml` | The input config plus the resolved λ, utilisations, A, ρ(A), α, scaled θ, seeds and package versions, and a `scenarios` set-up block for each scenario (λ, ρ(A), α, analytic and achieved utilisation, shock window, blocking events). |
| `r_change_vs_theta.png`, `bp_top2_gap_hist.png`, `contested_share_by_node.png` | The figures. |

### Per-epoch log columns

- `rep, time, node`: the replication, the simulation time in minutes, and the decision node.
- `triage_level`: r*, the most urgent level waiting. `n_feasible`: |F|. `n_distinct_next`: the number of distinct next nodes in F.
- `contested`: |F| ≥ 2 and at least 2 distinct next nodes.
- `bp_choice`: the class id chosen by PlainBP.
- `bp_top2_gap`: |P₁ − P₂| for the two best BP candidates (NaN when |F| = 1).
- `cent_gap_top2_k{κ}`: |ĉ_j1 − ĉ_j2| for the next nodes of those two candidates. There is one column per κ, because ĉ depends on κ.
- `q_total`, `q_<node>`: the aggregate waiting-queue vector q(t) and its total.
- `hop_gap_top2_k{κ}`: |ĥ_j1 − ĥ_j2| for the top-2 BP candidates.
- `oh_k{κ}_h{θ_h}`, `changed_oh_…`, `katzvs1h_…`: the OneHopBP choice, whether it differs from BP, and whether it differs from the KatzBP reference (κ = 0.5, θ = 0.1).
- `katz_k{κ}_{θ}`, `changed_k{κ}_{θ}`: the KatzBP choice and whether it differs from BP. Absolute θ columns are labelled `abs<θ>`; pressure-scaled θ columns are labelled `eta<η>`, and their values are listed in the report and in the resolved config.

## Metric definitions

Each metric is computed once per replication. The report shows the mean ± 95% CI across
replications (Student-t with n − 1 df).

- **Contested share**: contested epochs ÷ all decision epochs, overall and per node. Only contested epochs can ever flip.
- **R_change(κ, θ)**: epochs where KatzBP ≠ BP ÷ contested epochs, overall and per node. A change outside a contested epoch is impossible (the report checks that the count is 0).
- **bp_top2_gap / cent_gap_top2**: the median and IQR over contested epochs.
- **Dominance ratio(κ, θ)**: the median over contested epochs of θ·cent_gap_top2 ÷ bp_top2_gap. A value below 1 means the typical BP pressure difference outweighs the largest possible centrality shift, so BP still dominates.
- **BP tie share**: the share of contested epochs with bp_top2_gap = 0 exactly, overall and per node.
- **R_change excluding ties**: R_change restricted to contested epochs with bp_top2_gap > 0. It counts only real overrides of BP pressure, not tie-breaking, and it is the basis of the ACTIVE RANGE.
- **Candidate set**: the named (κ, θ) set from `katz.candidate_set` (κ = 0.5, θ ∈ {0.02, 0.1, 0.25}). The report and the console print its row from both R_change tables, together with its dominance ratio and status.
- **Utilisation**: the analytic value (λ × visits × E[S] ÷ servers) against the achieved post-warm-up busy fraction.
- **Decision time**: mean wall-clock µs per `KatzBP.select` call (Katz solve + argmax), with PlainBP as a reference.

### Katz vs one-hop

- **1H-BP** uses h = q + α Aᵀq, the first two terms of the Katz series, with the same α = κ/ρ(A). The bounded score is ĥ = h / (1 + max|h|), and W_c = μ_i^c (Q_i^c − Q_j^c) + θ_h (ĥ_i − ĥ_j). The tie-break is the same as for the other policies. It is shadow-evaluated at κ = 0.5 for θ_h ∈ {0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1} (config `one_hop:`).
- **R_change(1H vs BP)** and **R_change(Katz vs 1H)** count contested epochs where the two choices differ. The Katz side is the reference KatzBP at κ = 0.5, θ = 0.1. Both are reported including and excluding exact BP ties.
- **Dominance ratio (1H)**: the median over non-tied contested epochs of θ_h·|ĥ_j1 − ĥ_j2| ÷ bp_top2_gap.
- **Matched θ_h**: the θ_h in the grid whose mean 1H dominance ratio is closest, by absolute difference, to the reference KatzBP dominance ratio.
- **Verdict**: if Katz vs 1H excluding ties is below 5% at the matched θ_h → `KATZ ≈ ONE-HOP: full-network signal adds few distinct decisions.` Otherwise the report prints `KATZ ≠ ONE-HOP` with the disagreement rate.

### Scenario comparison

Each scenario runs the same replications and seeds under PlainBP, with FCFS at registration and triage. It shadow-evaluates one KatzBP (κ = 0.5, θ = 0.1) and one OneHopBP (κ = 0.5, θ_h = 0.1, the matched value). All settings are in `scenarios:` in the config; model overrides are deep-merged into `model:`, with dicts merged key by key and lists replaced whole.

| Scenario | Set-up |
|---|---|
| A baseline | The current config. |
| B high load | λ is recalibrated so the bottleneck has ρ = 0.95. |
| C lab shock | From t0 = warm-up + 2 days, lab services take twice as long (rate × 0.5) for 4 h. It is reported for the whole run and for the window [t0, t0 + 8 h] (C-window). |
| D blocking | Blocking-after-service, plus an `admission` node (13 beds, lognormal mean 240 min, CV 0.8, PROVISIONAL) that P5 and P6 visit after treatment. Waiting rooms hold at most lab 6, imaging 6, treatment 6, admission 4 (PROVISIONAL). A, ρ(A) and α are recomputed with admission, and the policy applies there too. |

For each row the report gives the contested share and three disagreement rates, including and excluding exact BP ties: Katz vs BP, 1H vs BP, and **Katz vs 1H** (the key metric), plus Katz vs 1H by node. Per-row verdict: `KATZ DISTINCT` if Katz vs 1H excluding ties ≥ 5%, else `KATZ ≈ ONE-HOP`. Scenario A shadows the same trajectory as the main pilot, so the report checks that its Katz vs 1H value equals the main pilot's value at θ_h = 0.1.

### Verdict rules (evaluated on the across-replication means)

1. Overall contested share < 2% → `MODEL: insufficient branching/competition — revise ED topology or class mix.`
2. R_change < 5% for every (κ, θ) → `BONUS TOO WEAK — rescale θ using the pressure gap.`
3. For each κ, the θ values with R_change **excluding ties** ≥ 5% and dominance ratio < 1 → `ACTIVE RANGE (excluding ties): θ ∈ [min, max], κ = …`
4. For each κ, the θ values with R_change including ties ≥ 5% but excluding ties < 5% → `TIE-BREAK ONLY`.
5. Any θ with R_change (including ties) > 50% → flagged as `centrality overrides BP — too strong.`

If none of rules 1–3 fires, the report prints `NO ACTIVE RANGE`. Rules 2 and 5 still use R_change including ties, as in the original spec. The thresholds are set in `verdict:` in the config.

## Tests

| Test | What it checks |
|---|---|
| `test_theta_zero_matches_plain_bp` | With θ = 0, KatzBP equals PlainBP on 1,000 random states for each κ. |
| `test_c_hat_bounded` | ĉ ∈ [0, 1) for random non-negative q, and ĉ(sink) = 0. |
| `test_katz_solve_matches_truncated_series` | The solve equals Σ αⁿ(Aᵀ)ⁿq to within 1e-8. |
| `test_zero_spectral_radius_fails_loudly` | An acyclic A raises an error. |
| `test_choice_is_most_urgent_level`, `test_triage_never_violated` | Triage order holds on random states and at every policy-node epoch in the simulation, under both BP and a strong KatzBP. FCFS nodes never produce decision epochs. |
| `test_fcfs_serves_earliest_join_regardless_of_level`, `test_fcfs_nodes_serve_in_arrival_order` | FCFS picks the earliest queue join, and registration and triage start service in join order. |
| `test_one_hop_equals_first_two_katz_terms` | h equals Σ_{n=0..1} αⁿ(Aᵀ)ⁿq, and ĥ ∈ [0, 1) with ĥ(sink) = 0. |
| `test_one_hop_theta_zero_matches_plain_bp_and_respects_triage` | With θ_h = 0, 1H-BP equals PlainBP, and it never breaks triage order. |
| `test_blocking_never_exceeds_buffer_and_conserves_patients` | In scenario D (config buffers, and buffers of 1), no waiting room ever exceeds its capacity, blocking does occur, and arrivals = departures + in system, under both BP and KatzBP. |
| `test_shock_changes_only_lab_service_within_window` | Only lab services starting in [t0, t0 + 4 h) are scaled (×2, from the same random draw). Everything else is unscaled, and the run is identical to the no-shock run before t0. |
| `test_scenario_d_adds_admission_to_network`, `test_deep_merge_…` | Admission is in A with capacity 4 and consultation stays the bottleneck; config overrides merge as documented. |
| `test_verdict.py` | The ACTIVE RANGE ignores tie-only changes, TIE-BREAK ONLY is labelled, `NO ACTIVE RANGE` fires when only ties change, and the matched θ_h selection and the KATZ ≈ / ≠ ONE-HOP rule behave as specified. |
| `test_patient_conservation` | arrivals = departures + in system. |
| `test_mm1_tandem_matches_jackson` | A 2-node M/M/1 tandem at ρ = 0.7 matches L = ρ/(1 − ρ) within 5%. The current result is 2.310 against 2.333. |
| `test_same_seed_same_summary` | The same seeds give an identical summary. Timing metrics are excluded because they are wall-clock. |
| Extra | Uncontested epochs never change; A is sub-stochastic; the calibration hits the target ρ. |

## Assumptions (simplest choice where the spec was open)

1. **Per-stage buffers.** P3 and P4 visit consultation twice, so a class (k, p) at consultation would have two different next nodes. Each (k, p, visit stage) queue is therefore a separate buffer, which is standard for multiclass networks. The policy's `class_id` is the buffer id, and Q_j^c is the queue of the same patient class's next-stage buffer at j.
2. **Queues count waiting patients only.** Q_i^c, Q_j^c and q_i exclude patients in service.
3. **Decision epoch.** An epoch is a server release at a policy node that finds patients waiting. Registration and triage are FCFS, so they have no epochs; their queues still count in q(t) for the Katz score. A patient who arrives at an idle server starts service immediately, which involves no choice because the work-conserving system has nobody else waiting. That dispatch is not logged as an epoch.
4. **Routing before deciding.** On a release, the finished patient joins the next node before the decision at the releasing node. The patient's move therefore shows in q(t) and in Q_j^c.
5. **The triage level is known from arrival.** The class is sampled at arrival. Registration and triage ignore it, because they are FCFS, and triage priority applies only at policy nodes.
6. **HOL time.** The tie-break uses the time the head-of-line patient joined that buffer's queue.
7. **Run length.** There are 24 h of warm-up plus 7 observed days, so each replication simulates 8 days in total.
8. **Class-specific service time.** The class mean is the node mean × `level_multiplier[k]`, × `revisit_multiplier` on a repeat visit (re-review). The CV is per node.
9. **Class mix** is given as Pr(k) × Pr(p | k) and printed as a joint table in `resolved_config.yaml`.
10. **Pressure-scaled θ.** A calibration pass of replication 1, using the same seed and trajectory, gives the median gap before the full runs. If that median were 0 because of ties, the median of the non-zero gaps would be used and the report would say so.
11. **Dominance ratio** excludes epochs with bp_top2_gap = 0, where the ratio is undefined. The tie share is reported separately.
12. **RNG streams.** Each stream is seeded by `SeedSequence([seed, stream_id])`, with fixed ids: arrivals 0, class assignment 1, policy 2, shadow 3, and service at node i 100 + i. This keeps common random numbers aligned when policies change later.
13. **The epoch log is CSV (gzip).** This avoids adding a parquet dependency.
14. **The active range uses R_change excluding ties** (≥ 5%) and a dominance ratio < 1. θ values above 50% R_change are flagged separately rather than removed from the range.
15. **FCFS order** is the time each patient joined that node's queue, with the class id breaking exact ties. Pure FCFS does not apply triage priority.
16. **Candidate θ in the grid.** A candidate θ that is missing from `theta_abs` (here 0.02) is added to the absolute grid for every κ, so it can be read from the same tables. In this config, θ = 0.02 also happens to equal the pressure-scaled η = 0.5 value.
17. **Matched θ_h** uses the absolute difference between mean dominance ratios. On an exact tie, the smaller θ_h wins.
18. **The KATZ ≠ ONE-HOP line** is an addition to the spec. It is printed when the ≈ condition fails, so the verdict always states the comparison.
19. **The KatzBP reference must already be in the Katz grid.** The pilot fails loudly if it is not, and the one-hop results are not included in the KatzBP decision-time metric.
20. **Admission beds (scenario D).** The spec's 6 beds would run admission at ρ = 1.72 at baseline λ. On your instruction, D uses 13 beds (ρ ≈ 0.79), so consultation stays the bottleneck and λ equals A's. The level multiplier applies to admission like every other node.
21. **Buffers** limit the number *waiting* at a node, excluding patients in service. A patient can always enter if a server is idle.
22. **Blocking.** A blocked patient keeps its server until a waiting place frees downstream. Blocked patients move in the order they were blocked, and the freed server's decision epoch happens only when the patient actually leaves. Blocked time counts as busy in achieved utilisation.
23. **Shock.** The shock applies to lab services that *start* inside [t0, t0 + 4 h): the same random draw is multiplied by 1/0.5. Services already in progress are not changed. Policies keep using the nominal μ, so the scheduler is not told about the shock.
24. **Shock window.** The C-window row uses the epochs with t0 ≤ time < t0 + 8 h. It has only about 36 contested epochs per replication, so its CIs are wide.
