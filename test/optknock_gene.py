# -*- coding: utf-8 -*-
"""
Gene-level OptKnock on iCG875v34  -- target product: L-tryptophan (EX_trp__L_e)

Core semantics (gene level != reaction level is the key):
  Reaction r is blocked  <=>  in the GPR's CNF there exists some OR clause C with C subset of the knockout set
  -> isoenzymes must be knocked out as a whole group; for a complex (AND), knocking out any subunit breaks it

Criterion (strict bilevel, excluding carbon-overflow false positives):
  mu*(S)     = max growth after knocking out S
  P_ko_max   = max Trp under fixed mu*      P_ko_min = min Trp under the same condition
  P_wt_max/min = upper/lower bound of Trp for the wild type at the same mu*
  OptKnock gain coupling : P_ko_max > P_wt_max + tol   (knockout forces the cell to overproduce Trp)
  RobustKnock forced coupling: P_ko_min > tol          (under optimal growth the cell has to produce Trp)
"""
import sys, time, itertools, re, warnings
import cobra, pandas as pd
warnings.filterwarnings("ignore")

def log(*a): print(*a, flush=True)      # key: flush, otherwise progress is invisible when redirecting

# ====================== tunable parameters ======================
MODEL     = "iCG875.xml"
OUT       = "supp_OptKnock_gene_Trp.csv"
GLC       = 0.467          # mmol/gDW/h, set by this project
TARGET    = 'EX_trp__L_e'
BIO       = 'CG_biomass_cgl_ATCC13032'
GLC_ID    = 'EX_glc_e'
MU_FRAC   = 0.50           # growth requirement = 50% of WT mu_max (strict tier; for a looser tier set 0.10)
K_MAX     = 3
HUB_N     = 25             # reactions count > this value -> treated as a hub metabolite, not traversed when tracing back
MAX_SETS  = 20000          # safety upper limit for enumeration
TOL       = 1e-5
# knocking these out is equivalent to switching off the Trp outlet/synthesis, not treated as candidate strategies
TRP_ESSENTIAL = {'TRPS1','TRPS2','TRPS3','TRPt2rpp','TRPabcpp','EX_trp__L_e'}
# ==============================================================

# ---------- GPR -> CNF (correct version) ----------
# AND: clause lists concatenated in parallel ; OR: clauses cross-merged
def _tok(s):
    return [t for t in s.replace('(', ' ( ').replace(')', ' ) ').split() if t]
def _atom(t, i):
    if i < len(t) and t[i] == '(':
        i += 1; sub, i = _and(t, i)
        if i < len(t) and t[i] == ')': i += 1
        return sub, i
    return [frozenset([t[i]])], i + 1
def _or(t, i):
    left, i = _atom(t, i)
    while i < len(t) and t[i].lower() == 'or':
        i += 1; right, i = _atom(t, i)
        left = [a | b for a in left for b in right]      # cross merge
    return left, i
def _and(t, i):
    left, i = _or(t, i)
    while i < len(t) and t[i].lower() == 'and':
        i += 1; right, i = _or(t, i)
        left = left + right                              # parallel concatenation
    return left, i
def gpr_to_cnf(rule):
    if not rule or not rule.strip(): return []
    try:
        cl, _ = _and(_tok(rule), 0); return cl
    except Exception:
        gs = re.findall(r'Cgl\d+', rule)
        if not gs: return []
        return [frozenset(gs)] if 'or' in rule.lower() else [frozenset([g]) for g in gs]

# ---------- self-check: the parser must pass ----------
_assert = lambda r, exp: (lambda g: None if [sorted(c) for c in g] == exp
                          else sys.exit(f"CNF parser self-check failed: {r} -> {g}"))(gpr_to_cnf(r))
_assert('Cgl3034 and Cgl3035',          [['Cgl3034'], ['Cgl3035']])
_assert('Cgl0700 or Cgl0831 or Cgl2870',[['Cgl0700','Cgl0831','Cgl2870']])
_assert('(Cgl1 and Cgl2) or Cgl3',      [['Cgl1','Cgl3'], ['Cgl2','Cgl3']])
log("[OK] GPR->CNF parser self-check passed")

m = cobra.io.read_sbml_model(MODEL)
log(f"model: {len(m.reactions)} reactions / {len(m.metabolites)} metabolites / {len(m.genes)} genes")
exch = set(r.id for r in m.exchanges)
m.reactions.get_by_id(GLC_ID).bounds = (-GLC, -GLC)
m.objective = BIO
mu_max_wt = m.optimize().objective_value
mu_min = MU_FRAC * mu_max_wt
log(f"glc={GLC} -> WT mu_max = {mu_max_wt:.5f} | growth threshold mu_min = {mu_min:.5f}")

# reference: the Trp upper bound under mu_min with no knockout (= carbon-overflow ceiling, used to judge whether a strategy is merely overflow)
ref = m.copy(); ref.reactions.get_by_id(BIO).lower_bound = mu_min
ref.objective = TARGET
overflow_ceiling = ref.optimize().objective_value
log(f"[reference] no knockout, max Trp under mu>=mu_min = {overflow_ceiling:.5f}  <- overflow ceiling, results should not merely approach it")

# ---------- candidate pool: Trp upstream (excluding hubs) + central metabolism ----------
HUB = {met.id for met in m.metabolites if len(met.reactions) > HUB_N}
def upstream(seeds, depth=3):
    seen, front = set(seeds), set(seeds)
    for _ in range(depth):
        nxt = set()
        for rid in front:
            for met, c in m.reactions.get_by_id(rid).metabolites.items():
                if c <= 0 or met.id in HUB or met.id in exch: continue
                for pr in met.reactions:
                    if pr.id in exch or pr.id in seen: continue
                    if pr.metabolites.get(met, 0) > 0: nxt.add(pr.id)
        nxt -= seen
        if not nxt: break
        seen |= nxt; front = nxt
    return seen
rt = m.reactions.get_by_id(TARGET)
seeds = [rt.id] + [pr.id for met in rt.metabolites for pr in met.reactions if pr.id not in exch]
path = upstream(sorted(set(seeds)), 3)
CENTRAL = ['PFK','PYK','GAPD','PGK','PGM','ENO','FBA','TPI','ZWF','GND','PGL','RPE','RPI','TKT','TALA',
 'CS','ACONT','IDH','AKGDH','SUCOAS','SDH','FUM','MDH','PPC','PPCK','PDH','ME2','G6PDH2r','HEX',
 'HISTD','HISTP','HISTF','HISTA','ATPPRT','PRPPS','IGPS','TRPS1','TRPA','CHORS','CHORM','SHKC','SHSK',
 'DHNPA','ASPTA','GLUD','GLNS','SERAT','LYSD','DAPDC','AROD','AROB','AROC','TYRTA','PHEATA','SHMT','GLYAT']
central = {r.id for r in m.reactions if r.id not in exch
           and any(k.upper() in r.id.upper() for k in CENTRAL)}
cand_rxns = sorted((path | central) - exch)
CNF = {rid: gpr_to_cnf(m.reactions.get_by_id(rid).gene_reaction_rule) for rid in cand_rxns}
genes = sorted({g for rid in cand_rxns for cl in CNF[rid] for g in cl})
log(f"candidates: {len(cand_rxns)} reactions / {len(genes)} genes  (Trp upstream {len(path)} + central metabolism {len(central)})")

# ---------- clause-driven enumeration: only generate <=K gene sets that 'indeed block some reaction' ----------
sets = set()
for rid in cand_rxns:
    for cl in CNF[rid]:
        if 1 <= len(cl) <= K_MAX: sets.add(frozenset(cl))
base = sorted(sets, key=lambda s: (len(s), sorted(s)))
for a, b in itertools.combinations(base, 2):
    u = a | b
    if 2 <= len(u) <= K_MAX: sets.add(frozenset(u))
cand_sets = sorted(sets, key=lambda s: (len(s), sorted(s)))
log(f"candidate knockout sets: {len(cand_sets)}  (K=1: {sum(1 for s in cand_sets if len(s)==1)}, "
    f"K=2: {sum(1 for s in cand_sets if len(s)==2)}, K=3: {sum(1 for s in cand_sets if len(s)==3)})")
if len(cand_sets) > MAX_SETS:
    import random; random.seed(0)
    log(f"[WARN] exceeds upper limit {MAX_SETS}, sampling randomly"); cand_sets = random.sample(cand_sets, MAX_SETS)

# ---------- bilevel validation ----------
wt_cache = {}
def wt_bounds(mu):
    key = round(mu, 6)
    if key in wt_cache: return wt_cache[key]
    w = m.copy(); w.reactions.get_by_id(BIO).bounds = (mu - TOL, mu + TOL)
    w.objective = TARGET; w.solver.objective.direction = 'max'
    mx = w.optimize().objective_value
    w.solver.objective.direction = 'min'
    mn = w.optimize().objective_value
    wt_cache[key] = (mx, mn); return wt_cache[key]

stat = dict(no_block=0, kills_trp=0, no_growth=0, ok=0)
rows, t0 = [], time.time()
for idx, ks in enumerate(cand_sets, 1):
    blocked = [rid for rid in cand_rxns if CNF[rid] and any(cl <= ks for cl in CNF[rid])]
    if not blocked:                       stat['no_block']  += 1; continue
    if set(blocked) & TRP_ESSENTIAL:      stat['kills_trp'] += 1; continue
    mp = m.copy()
    for rid in blocked: mp.reactions.get_by_id(rid).bounds = (0, 0)
    mp.objective = BIO
    s = mp.optimize()
    if s.status != 'optimal' or s.objective_value < mu_min - 1e-6:
        stat['no_growth'] += 1; continue
    mu = s.objective_value
    mp.reactions.get_by_id(BIO).bounds = (mu - TOL, mu + TOL)
    mp.objective = TARGET; mp.solver.objective.direction = 'max'
    pk_max = mp.optimize().objective_value
    mp.solver.objective.direction = 'min'
    pk_min = mp.optimize().objective_value
    pw_max, pw_min = wt_bounds(mu)
    stat['ok'] += 1
    gain = pk_max - pw_max
    rows.append(dict(gene_knockouts='|'.join(sorted(ks)), n_gene_ko=len(ks),
                     blocked_reactions='|'.join(blocked), n_blocked=len(blocked),
                     mu_star=round(mu, 6), mu_ratio=round(mu / mu_max_wt, 4),
                     p_trp_ko_max=round(pk_max, 6), p_trp_ko_min=round(pk_min, 6),
                     p_trp_wt_max=round(pw_max, 6), p_trp_wt_min=round(pw_min, 6),
                     gain_vs_wt=round(gain, 6),
                     yield_g_per_glc=round(pk_max / GLC, 4),
                     optknock_gain_coupled=bool(gain > 1e-4),
                     robust_growth_coupled=bool(pk_min > 1e-6)))
    if stat['ok'] > 0 and stat['ok'] % 200 == 0:
        log(f"  progress {idx}/{len(cand_sets)} | valid {stat['ok']} | coupled found "
            f"{sum(r['optknock_gain_coupled'] for r in rows)} | {time.time()-t0:.0f}s")

df = pd.DataFrame(rows)
log(f"\nstats: valid={stat['ok']}  block no reaction={stat['no_block']}  "
    f"hurts Trp pathway={stat['kills_trp']}  growth below threshold={stat['no_growth']}  elapsed {time.time()-t0:.0f}s")

if len(df):
    df = df.sort_values(['robust_growth_coupled','optknock_gain_coupled','gain_vs_wt'],
                        ascending=False)
    df.to_csv(OUT, index=False)
    log(f"[OK] saved {OUT} ({len(df)} rows)")
    log(f"  OptKnock gain-coupled strategies: {int(df.optknock_gain_coupled.sum())}")
    log(f"  RobustKnock forced-coupled strategies: {int(df.robust_growth_coupled.sum())}")
    log("\n=== top 15 candidates with strongest coupling criterion / largest gain ===")
    cols = ['gene_knockouts','n_gene_ko','mu_star','p_trp_ko_max','p_trp_wt_max','gain_vs_wt',
            'yield_g_per_glc','optknock_gain_coupled','robust_growth_coupled']
    log(df[cols].head(15).to_string(index=False))
    log("\n=== suggested experimental priority (gain>0 and growth retention>50% mu_max) ===")
    pri = df[(df.gain_vs_wt > 1e-4) & (df.mu_ratio >= 0.5)]
    log(pri[cols].to_string(index=False) if len(pri) else "  (none - see interpretation below)")
else:
    pd.DataFrame(columns=['gene_knockouts']).to_csv(OUT, index=False)
    log("[WARN] no valid knockout set produced")
