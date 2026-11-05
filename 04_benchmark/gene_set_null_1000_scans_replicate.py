#!/usr/bin/env python
"""gene_set_null_1000_scans_replicate.py (2026-10-19, package geneset_null_1000_20261019; revised the same day after two reviews)

One replicate of the STAAR stage of the permuted gene-set null on the stratified design (covset strat_project): a provenance
wrapper around the package copy of rare_variant_tools_smmat_genesis_staar.R (--which staar --covset strat_project), mirroring the status-JSON
convention of five_analysis_execution_20261016/external_tools/corrected_v2/run_tool_task.py (one status JSON per task with driver and
input sha256s, native phenotype source and sha256, seed, accepted flag; atomic writes; no silent overwrite; accepted results are
never recomputed).

Steps
  1. <status-dir>/<key>_status.json (nothing is written before these checks pass):
       accepted      -> result sha256s re-checked, "Already finished", exit 0
       needs_audit   -> exit 2, nothing touched (after the audit a human moves the status and attempt files to an audited/ folder)
       running       -> LOCK: if the recorded SLURM job is still queued or running (squeue -j <id>) the replicate is in flight
                        and this task exits 3 without touching anything; when squeue cannot be asked, a 'running' status younger
                        than the time limit (+15 min) is also a lock; otherwise the earlier task died without a handler (node
                        failure, SIGKILL) and its status is moved aside as <key>_status.<stamp>.running_dead.json
       failed / timeout / terminated -> moved aside as <key>_status.<stamp>.<status>.json
     Every moved-aside file is one prior attempt; with two or more prior attempts the task refuses to run (needs_audit, exit 2),
     so repeated time-outs are counted and become visible instead of being resubmitted forever.
  2. Partial outputs <key>_genes.tsv.gz.part.<jobid> / <key>_meta.json.part.<jobid> (the R copy writes them with the job id of
     the writing task): if their job is still listed by squeue the replicate is in flight -> exit 3, nothing touched; partial
     files of dead jobs are removed as stale (recorded in the status JSON).
  3. Label file ystar_rep<NNN>.csv: regenerated in memory from the bignull npz (make_ystar_labels.label_frame) and compared byte
     for byte with the file on disk (written atomically if absent). The npz sha256 and seed are recorded.
  4. Outputs <out-dir>/strat_project_rep<NNN>_{genes.tsv.gz,meta.json}: if both exist they are validated and, when valid, accepted
     without recomputation (r_stage "existing"); invalid outputs are moved aside (*.rejected_<stamp>) and recomputed.
  5. Rscript rare_variant_tools_smmat_genesis_staar.R ... --which staar; its stdout+stderr are streamed to <status-dir>/<key>_R.log (appended per
     attempt) and to the task log. Acceptance requires ALL of: exit code 0; the literal line 'STAAR fit_null_glm: converged TRUE'
     in the R log (the R copy also stops on a non-converged null) AND staar.converged == true in meta.json; meta.json
     slurm_job_id == this task's SLURM_JOB_ID (ownership: the final files were renamed into place by THIS job); validation
     against the frozen reference replicate strat_project_rep001 (20,786 analyzed, 7,598 cases, design rank 85, mask coding,
     26,768 genes in the same order with the same n_var and cMAC; the SAME set of genes with a finite STAAR-O, 21,853, with
     ZERO tolerance: a gene finite in the reference but NA here is a per-gene STAAR error swallowed by the R tryCatch and
     rejects the replicate, the gene names being listed in the status JSON; every p_staar_* column finite exactly where
     STAAR-O is; p in [0, 1]).
  6. Signals: SIGUSR1 (SLURM sends it to the batch shell 600 s before the time limit, --signal=B:USR1@600, and the sbatch script
     forwards it) is recorded in the status JSON; SIGTERM/SIGINT raise an exception so that the except/finally blocks run: the
     R child is terminated and the status is written as 'timeout' (a USR1 warning was received, or the elapsed time is within
     15 min of the time limit) or 'terminated' (scancel, node drain) before SLURM's SIGKILL (KillWait 30 s).
  7. --compare-existing: the produced p_staar_* columns are compared with the frozen result file of the same replicate (the
     validation gate of the package: same NaN pattern and max |delta log10 p| <= 1e-4); report in <status-dir>/<key>_compare.json.

Exit codes: 0 accepted, 1 failed, 2 needs_audit, 3 locked / in flight (nothing done), 124 timeout, 143 terminated.

Usage
  gene_set_null_1000_scans_replicate.py --rep 101 --time-limit 10:00:00              (array task; out-dir = outputs/plasmode/rv_tools/genome/tools)
  gene_set_null_1000_scans_replicate.py --rep 1 --out-dir <pkg>/validation --status-dir <pkg>/validation --compare-existing   (validation gate)
  gene_set_null_1000_scans_replicate.py --rep 7 --validate-only                      (read-only validation report of an existing replicate; nothing is run)
  gene_set_null_1000_scans_replicate.py --rep 101 --out-dir <tmp> --status-dir <tmp> --ydir <tmp> --dry-r 60   (TEST HOOK: a sleeping stand-in for Rscript;
                                                          refused in the frozen directories; used only to exercise the signal path)
"""
import argparse, hashlib, json, os, re, resource, signal, socket, subprocess, sys, time, traceback
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_ystar_labels as myl  # noqa: E402

PACKAGE = "geneset_null_1000_20261019"
ROOT = Path("${ADSP_ROOT}")
G = ROOT / "outputs/plasmode/rv_tools/genome"
TOOLS = G / "tools"
YDIR = G / "ystar"
GDS = G / "coding_genome.seq.gds"
GROUPS = G / "groups/gene_variants.csv.gz"
POOL = ROOT / "outputs/plasmode/plasmode_pool.csv"
REF_STEM = TOOLS / "strat_project_rep001"          # frozen reference replicate (2026-10-02 run)
RSCRIPT = HERE / "rare_variant_tools_smmat_genesis_staar.R"
COVSET, WHICH, TOOL = "strat_project", "staar", "staar"
PCOLS = ["p_staar_O", "p_staar_S_1_25", "p_staar_B_1_25", "p_staar_A_1_25", "p_staar_S_1_1", "p_staar_B_1_1", "p_staar_A_1_1"]
META_KEYS = ["covset", "n_analyzed", "n_cases", "n_design_cols", "design_rank", "dropped_columns", "mask"]
DRIVER_FILES = ["gene_set_null_1000_scans_replicate.py", "make_ystar_labels.py", "rare_variant_tools_smmat_genesis_staar.R", "geneset_null_staar.sbatch", "geneset_null_validate.sbatch",
                "geneset_null_ystar.sbatch", "geneset_null_summary.sbatch", "gene_set_null_1000_scans_summary.py", "completion_report.py"]
LOG10_TOL = 1e-4
DEFAULT_TIME_LIMIT = "10:00:00"
TIMEOUT_MARGIN = 900                     # s: a SIGTERM this close to the time limit is a time-out even without the USR1 warning
ATTEMPT_LABELS = ("failed", "timeout", "terminated", "running_dead")
CONVERGED_RE = re.compile(r"STAAR fit_null_glm: converged (TRUE|FALSE)")
STATE = {"usr1_at": None, "signal": None, "child": None}


class Terminated(Exception):
    pass


def sha(path):
    return myl.sha256(path)


def write_json(path, obj, tmp_suffix=".pending.json"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(tmp_suffix)
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, path)


def stamp():
    return time.strftime("%Y%m%dT%H%M%S")


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def outputs(stem):
    return Path(str(stem) + "_genes.tsv.gz"), Path(str(stem) + "_meta.json")


def parse_slurm_time(s):
    """squeue %l / sbatch --time formats ('10:00:00', '1-00:00:00', '1-12', '30:00' = mm:ss, '45' = minutes) -> seconds; None if unlimited/unparsable."""
    s = (s or "").strip()
    if not s or s.upper() in ("UNLIMITED", "NOT_SET", "N/A", "INVALID"):
        return None
    try:
        days = 0
        if "-" in s:
            d, s = s.split("-", 1)
            days = int(d)
        parts = [int(x) for x in s.split(":")]
        if len(parts) == 3:
            h, m, sec = parts
        elif len(parts) == 2:
            h, m, sec = (parts[0], parts[1], 0) if days else (0, parts[0], parts[1])
        elif len(parts) == 1:
            h, m, sec = (parts[0], 0, 0) if days else (0, parts[0], 0)
        else:
            return None
        return days * 86400 + h * 3600 + m * 60 + sec
    except ValueError:
        return None


def job_alive(jobid):
    """Is SLURM job `jobid` still queued, running or completing? True / False, or None when squeue cannot be asked."""
    if not jobid:
        return None
    try:
        r = subprocess.run(["squeue", "-h", "-j", str(jobid), "-o", "%T"], capture_output=True, text=True, timeout=120)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return None
    if r.returncode == 0:
        return bool(r.stdout.strip())
    if "Invalid job id" in (r.stderr or ""):
        return False
    return None


def part_files(gf, mf):
    return sorted(set(gf.parent.glob(gf.name + ".part*")) | set(mf.parent.glob(mf.name + ".part*")))


def part_owner(p):
    m = re.search(r"\.part\.(\d+)$", p.name)
    return m.group(1) if m else None


def in_flight(owner, age, limit_s, my_job):
    """A lock holder is in flight when squeue lists its job, or (squeue unavailable) when it is younger than the time limit."""
    if owner and my_job and str(owner) == str(my_job):
        return False
    alive = job_alive(owner)
    if alive is True:
        return True
    if alive is None and age < limit_s + TIMEOUT_MARGIN:
        return True
    return False


def load_reference():
    gf, mf = outputs(REF_STEM)
    ref_meta = json.loads(mf.read_text())
    ref = pd.read_csv(gf, sep="\t", usecols=["gene", "chr", "n_var", "cMAC", "pos_min", "pos_max", "p_staar_O"])
    ref["finite_O"] = np.isfinite(ref.p_staar_O.to_numpy(float))
    return ref_meta, ref


def validate_outputs(stem, ref_meta, ref_genes, require_converged=True, expect_job_id=None):
    """Validate one replicate's outputs against the frozen reference. Returns (ok, report). Zero tolerance on the finite pattern."""
    gf, mf = outputs(stem)
    rep = dict(genes_file=str(gf), meta_file=str(mf), problems=[], warnings=[])
    if not (gf.exists() and mf.exists()):
        rep["problems"].append("missing output file(s)")
        return False, rep
    try:
        m = json.loads(mf.read_text())
    except Exception as e:
        rep["problems"].append(f"meta.json unreadable: {e}")
        return False, rep
    for k in META_KEYS:
        if m.get(k) != ref_meta.get(k):
            rep["problems"].append(f"{k}: {m.get(k)!r} != reference {ref_meta.get(k)!r}")
    st = m.get("staar") or {}
    n_ref = int(ref_meta["staar"]["n_genes"])
    rep.update(meta_n_analyzed=m.get("n_analyzed"), meta_n_cases=m.get("n_cases"), meta_design_rank=m.get("design_rank"),
               staar_n_genes=st.get("n_genes"), staar_seconds=st.get("seconds"), staar_lambda=st.get("lambda"), total_seconds=m.get("seconds"),
               staar_converged_meta=st.get("converged"), meta_slurm_job_id=m.get("slurm_job_id"), meta_package=m.get("package"))
    if "staar" not in m:
        rep["problems"].append("meta.json has no 'staar' block (STAAR stage did not run)")
    conv = st.get("converged")
    if conv is not True:
        if require_converged or conv is False:
            rep["problems"].append(f"meta staar.converged is {conv!r} (true required: the STAAR logistic null must have converged)")
        else:
            rep["warnings"].append("meta.json has no staar.converged key (frozen three-tool run of 2026-10-02; its log shows converged TRUE)")
    if expect_job_id is not None and str(m.get("slurm_job_id", "")) != str(expect_job_id):
        rep["problems"].append(f"meta slurm_job_id {m.get('slurm_job_id')!r} != this task's SLURM_JOB_ID {expect_job_id!r}: the final files were written by another job")
    try:
        d = pd.read_csv(gf, sep="\t")
    except Exception as e:
        rep["problems"].append(f"genes file unreadable: {e}")
        return False, rep
    rep["n_rows"] = int(d.shape[0])
    same_order = d.shape[0] == ref_genes.shape[0] and np.array_equal(d.gene.to_numpy(), ref_genes.gene.to_numpy())
    if not same_order:
        rep["problems"].append(f"gene list/order differs from the reference ({d.shape[0]} vs {ref_genes.shape[0]} rows)")
    else:
        for c in ("chr", "n_var", "cMAC", "pos_min", "pos_max"):
            if not np.array_equal(d[c].to_numpy(), ref_genes[c].to_numpy()):
                rep["problems"].append(f"column {c} differs from the reference")
    missing = [c for c in PCOLS + ["n_var_staar"] if c not in d.columns]
    if missing:
        rep["problems"].append(f"missing columns {missing}")
    else:
        P = d[PCOLS].to_numpy(float)
        fin = np.isfinite(P)
        fin_O = fin[:, 0]
        nfin = int(fin_O.sum())
        rep["n_finite_staar_O"] = nfin
        rep["n_finite_reference"] = n_ref
        if st.get("n_genes") is not None and int(st["n_genes"]) != nfin:
            rep["problems"].append(f"meta staar.n_genes {st['n_genes']} != finite STAAR-O count {nfin}")
        if same_order:
            fref = ref_genes.finite_O.to_numpy()
            na_only = ref_genes.gene.to_numpy()[fref & ~fin_O].tolist()      # finite in the reference, NA here: per-gene STAAR error
            fin_only = ref_genes.gene.to_numpy()[~fref & fin_O].tolist()
            rep["n_genes_na_here_finite_in_reference"] = len(na_only)
            rep["n_genes_finite_here_na_in_reference"] = len(fin_only)
            if na_only or fin_only:
                rep["genes_na_here_finite_in_reference"] = na_only[:50]
                rep["genes_finite_here_na_in_reference"] = fin_only[:50]
                rep["problems"].append(f"finite STAAR-O pattern differs from the reference (zero tolerance): {len(na_only)} gene(s) NA here but finite "
                                       f"in the reference (per-gene STAAR error swallowed by tryCatch in the R worker), {len(fin_only)} finite here but NA there")
        elif nfin != n_ref:
            rep["problems"].append(f"finite STAAR-O count {nfin} != reference {n_ref}")
        if not np.array_equal(fin.all(axis=1), fin_O):
            rep["problems"].append("NaN pattern differs between the p_staar_* columns (every column must be finite exactly where STAAR-O is)")
        if fin.any() and (np.any(P[fin] < 0) or np.any(P[fin] > 1)):
            rep["problems"].append("p values outside [0, 1]")
    return (not rep["problems"]), rep


def compare_with_frozen(stem, rep_no):
    """Compare the p_staar_* columns of <stem> with the frozen result of the same replicate in TOOLS."""
    new = pd.read_csv(outputs(stem)[0], sep="\t")
    old = pd.read_csv(outputs(TOOLS / f"strat_project_rep{rep_no:03d}")[0], sep="\t")
    out = dict(replicate=rep_no, new=str(outputs(stem)[0]), frozen=str(outputs(TOOLS / f"strat_project_rep{rep_no:03d}")[0]),
               same_gene_order=bool(np.array_equal(new.gene.to_numpy(), old.gene.to_numpy())), columns={})
    worst = 0.0
    for c in PCOLS + ["n_var_staar"]:
        a, b = new[c].to_numpy(float), old[c].to_numpy(float)
        fa, fb = np.isfinite(a), np.isfinite(b)
        same_nan = bool(np.array_equal(fa, fb))
        both = fa & fb
        if c == "n_var_staar":
            d = np.abs(a[both] - b[both])
        else:
            d = np.abs(np.log10(np.clip(a[both], 1e-300, 1)) - np.log10(np.clip(b[both], 1e-300, 1)))
        out["columns"][c] = {"n_compared": int(both.sum()), "same_nan_pattern": same_nan, "n_exact_equal": int(np.sum(a[both] == b[both])),
                             "max_abs_delta": float(d.max()) if d.size else 0.0, "q999_abs_delta": float(np.percentile(d, 99.9)) if d.size else 0.0,
                             "n_delta_gt_1e-8": int(np.sum(d > 1e-8))}
        if not same_nan:
            out["columns"][c]["n_nan_mismatch"] = int(np.sum(fa != fb))
        worst = max(worst, out["columns"][c]["max_abs_delta"] if c != "n_var_staar" else (1.0 if d.size and d.max() > 0 else 0.0))
    out["max_abs_delta_log10p"] = worst
    out["tolerance_log10p"] = LOG10_TOL
    out["pass"] = bool(out["same_gene_order"] and all(v["same_nan_pattern"] for v in out["columns"].values()) and worst <= LOG10_TOL)
    return out


def install_signals(meta, status):
    """USR1 = SLURM's pre-limit warning (recorded); TERM/INT = raise Terminated once so that except/finally run; repeats are ignored."""
    def on_usr1(signum, frame):
        STATE["usr1_at"] = time.time()
        meta.setdefault("signals", []).append(dict(signal="SIGUSR1", at=now(), note="SLURM time-limit warning (600 s before the limit)"))
        meta["time_limit_warning_at"] = now()
        try:
            write_json(status, meta, tmp_suffix=".pending_usr1.json")
        except Exception:
            pass
        sys.stdout.write(f"[run_rep] SIGUSR1 at {now()}: within 600 s of the time limit\n"); sys.stdout.flush()

    def on_term(signum, frame):
        if STATE["signal"] is not None:            # second TERM (SLURM + forwarded by the sbatch shell): already unwinding
            return
        STATE["signal"] = int(signum)
        meta.setdefault("signals", []).append(dict(signal=signal.Signals(signum).name, at=now()))
        raise Terminated(f"signal {signal.Signals(signum).name}")

    signal.signal(signal.SIGUSR1, on_usr1)
    signal.signal(signal.SIGTERM, on_term)
    signal.signal(signal.SIGINT, on_term)


def run_r(cmd, logfile):
    """Run the R driver, streaming stdout+stderr to `logfile` and the task log; returns (exit code, converged flag from the log)."""
    conv = None
    with open(logfile, "a") as lf:
        lf.write(f"==== {now()} job {os.getenv('SLURM_JOB_ID')} host {socket.gethostname()}: {' '.join(cmd)}\n"); lf.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        STATE["child"] = proc
        try:
            for line in proc.stdout:
                lf.write(line); lf.flush()
                sys.stdout.write(line); sys.stdout.flush()
                mm = CONVERGED_RE.search(line)
                if mm:
                    conv = (mm.group(1) == "TRUE")
            rc = proc.wait()
            lf.write(f"==== {now()} Rscript exit code {rc}\n")
        except Terminated:
            lf.write(f"==== {now()} wrapper received signal {STATE['signal']}: terminating Rscript\n"); lf.flush()
            try:
                proc.terminate(); proc.wait(timeout=15)
            except Exception:
                try:
                    proc.kill()
                except Exception:
                    pass
            raise
        finally:
            STATE["child"] = None
    return rc, conv


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--out-dir", default=str(TOOLS))
    ap.add_argument("--status-dir", default=str(HERE / "status"))
    ap.add_argument("--ydir", default=str(YDIR), help="label directory (the benchmark's by default; only the test hook changes it)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--time-limit", default=DEFAULT_TIME_LIMIT, help="this task's SLURM time limit (squeue %%l format): classifies a SIGTERM as a time-out and bounds the age of a stale lock")
    ap.add_argument("--compare-existing", action="store_true", help="validation gate: compare with the frozen result of the same replicate")
    ap.add_argument("--validate-only", action="store_true", help="read-only validation report of existing outputs; runs nothing")
    ap.add_argument("--dry-r", type=float, default=0, help="TEST HOOK: stand-in for Rscript that prints the converged line and sleeps this many seconds")
    a = ap.parse_args()
    assert 1 <= a.rep <= 1000
    out_dir, status_dir, ydir = Path(a.out_dir), Path(a.status_dir), Path(a.ydir)
    key = f"strat_project_rep{a.rep:03d}"
    stem = out_dir / key
    gf, mf = outputs(stem)
    in_frozen_dir = out_dir.resolve() == TOOLS.resolve()
    limit_s = parse_slurm_time(a.time_limit) or parse_slurm_time(DEFAULT_TIME_LIMIT)
    my_job = os.getenv("SLURM_JOB_ID")
    ref_meta, ref_genes = load_reference()
    if a.validate_only:
        ok, rep = validate_outputs(stem, ref_meta, ref_genes, require_converged=not (in_frozen_dir and a.rep <= 100))
        print(json.dumps(dict(key=key, valid=ok, **rep), indent=2, default=str))
        return 0 if ok else 1
    if in_frozen_dir and a.rep <= 100:
        raise SystemExit(f"replicates 1-100 are the frozen 2026-10-02 results in {TOOLS}; refusing to touch {key}")
    if a.compare_existing and in_frozen_dir:
        raise SystemExit("--compare-existing needs an --out-dir other than the frozen results directory")
    if a.dry_r and (in_frozen_dir or ydir.resolve() == YDIR.resolve() or status_dir.resolve() == (HERE / "status").resolve()):
        raise SystemExit("--dry-r is a test hook: it needs test --out-dir, --status-dir and --ydir outside the package and benchmark directories")
    status = status_dir / f"{key}_status.json"
    status_dir.mkdir(parents=True, exist_ok=True)
    # ---- 1. status checks: nothing is written before they pass
    if status.exists():
        old = json.loads(status.read_text())
        ost = old.get("status")
        if ost == "accepted":
            for f, h in ((old.get("result_genes"), old.get("result_genes_sha256")), (old.get("result_meta"), old.get("result_meta_sha256"))):
                assert f and Path(f).exists() and sha(f) == h, f"accepted status but result changed or vanished: {f}"
            print("Already finished", key, "accepted", flush=True)
            return 0
        if ost == "needs_audit":
            print(f"STATUS {key} needs_audit (unchanged): {old.get('reason')}; nothing run", flush=True)
            return 2
        label = ost or "unknown"
        if ost == "running":
            age = time.time() - float(old.get("start", 0) or 0)
            if in_flight(old.get("slurm_job_id"), age, limit_s, my_job):
                print(f"LOCKED: {key} is in flight in SLURM job {old.get('slurm_job_id')} (array {old.get('slurm_array_job_id')}_{old.get('slurm_array_task_id')}, "
                      f"host {old.get('hostname')}, started {old.get('started_at')}, {age / 3600:.2f} h ago); refusing to run a duplicate; nothing touched", flush=True)
                return 3
            label = "running_dead"
        moved = status.with_name(f"{key}_status.{stamp()}.{label}.json")
        os.replace(status, moved)
        print(f"prior status {ost} moved to {moved}", flush=True)
    prior = sorted(p for p in status_dir.glob(f"{key}_status.*.json") if p.name.split(".")[-2] in ATTEMPT_LABELS)
    if len(prior) >= 2:
        write_json(status, dict(key=key, package=PACKAGE, replicate=a.rep, status="needs_audit", accepted=False,
                                reason=f"{len(prior)} prior attempts ({', '.join(p.name.split('.')[-2] for p in prior)}); read the attempt files and the R log before any rerun",
                                prior_attempts=[str(p) for p in prior], r_log=str(status_dir / f"{key}_R.log"), finished=time.time(), finished_at=now()))
        print("STATUS", key, "needs_audit", flush=True)
        return 2
    # ---- 2. partial outputs of another task that is still in flight?
    parts = part_files(gf, mf)
    for p in parts:
        if in_flight(part_owner(p), time.time() - p.stat().st_mtime, limit_s, my_job):
            print(f"LOCKED: partial output {p} belongs to SLURM job {part_owner(p) or '(unknown; younger than the time limit)'}, still in flight; refusing to run a duplicate; nothing touched", flush=True)
            return 3
    meta = dict(key=key, package=PACKAGE, tool=TOOL, design=COVSET, covset=COVSET, which=WHICH, replicate=a.rep, status="running", accepted=False,
                attempt=len(prior) + 1, prior_attempts=[str(p) for p in prior], out_dir=str(out_dir), hostname=socket.gethostname(),
                slurm_job_id=my_job, slurm_array_job_id=os.getenv("SLURM_ARRAY_JOB_ID"), slurm_array_task_id=os.getenv("SLURM_ARRAY_TASK_ID"),
                time_limit=a.time_limit, time_limit_seconds=limit_s, start=time.time(), started_at=now(),
                driver_sha256={f: sha(HERE / f) for f in DRIVER_FILES if (HERE / f).exists()},
                input_sha256={"plasmode_pool.csv": sha(POOL), "coding_genome.seq.gds": sha(GDS), "gene_variants.csv.gz": sha(GROUPS)},
                reference_meta_sha256=sha(outputs(REF_STEM)[1]), reference_values={k: ref_meta.get(k) for k in META_KEYS} | {"staar_n_genes": ref_meta["staar"]["n_genes"]},
                r_log=str(status_dir / f"{key}_R.log"))
    write_json(status, meta)
    install_signals(meta, status)
    rc_final = 1
    try:
        # ---- 3. phenotype: regenerate in memory, compare byte for byte with the label file (write it if absent)
        pool = myl.pool_iid()
        df, info = myl.label_frame(a.rep, pool)
        ycsv = ydir / f"ystar_rep{a.rep:03d}.csv"
        want = myl.frame_bytes(df)
        if ycsv.exists():
            if ycsv.read_bytes() != want:
                raise RuntimeError(f"label file {ycsv} differs from the regenerated content of {info['npz']}")
            action = "existing_identical"
        else:
            ydir.mkdir(parents=True, exist_ok=True)
            myl.write_frame(df, ycsv)
            action = "written"
        meta.update(native_source=info["npz"], native_sha256=sha(info["npz"]), seed=info["seed"], n_cases_label=info["n_cases"],
                    phenotype=str(ycsv), phenotype_sha256=sha(ycsv), phenotype_action=action)
        write_json(status, meta)
        # ---- 4. existing outputs: accept if valid, else move aside and recompute
        r_stage = "run"
        if gf.exists() and mf.exists():
            ok, rep = validate_outputs(stem, ref_meta, ref_genes, require_converged=True)
            rlog = status_dir / f"{key}_R.log"
            if ok and rlog.exists() and "converged FALSE" in rlog.read_text():
                ok = False
                rep["problems"].append(f"{rlog} records 'converged FALSE'")
            if ok:
                r_stage = "existing"
                meta["validation"] = rep
                print("outputs exist and validate; accepted without recomputation", flush=True)
            else:
                st = stamp()
                for f in (gf, mf):
                    os.replace(f, f.with_name(f.name + f".rejected_{st}"))
                meta.setdefault("rejected", []).append(dict(stamp=st, problems=rep["problems"]))
                print(f"existing outputs rejected ({rep['problems']}); moved aside with suffix .rejected_{st}", flush=True)
        if r_stage == "run":
            for p in part_files(gf, mf):                 # dead jobs' partial files (the in-flight ones made us exit 3 above)
                p.unlink()
                meta.setdefault("removed_stale_parts", []).append(str(p))
                print("removed stale partial output of a dead job:", p, flush=True)
            if gf.exists() or mf.exists():
                raise RuntimeError(f"partial outputs present for {key} (genes {gf.exists()}, meta {mf.exists()}); audit before rerun")
            out_dir.mkdir(parents=True, exist_ok=True)
            cmd = ["Rscript", str(RSCRIPT), "--ystar", str(ycsv), "--covset", COVSET, "--gds", str(GDS), "--groups", str(GROUPS),
                   "--out", str(stem), "--which", WHICH, "--threads", str(a.threads)]
            if a.dry_r:
                cmd = [sys.executable, "-c", f"import sys, time; print('[test] STAAR fit_null_glm: converged TRUE (0 s)'); sys.stdout.flush(); time.sleep({a.dry_r}); print('[test] done')"]
            meta["r_command"] = " ".join(cmd)
            print("RUN", meta["r_command"], flush=True)
            t0 = time.time()
            rc, conv = run_r(cmd, status_dir / f"{key}_R.log")
            ru = resource.getrusage(resource.RUSAGE_CHILDREN)
            meta.update(r_exit_code=rc, r_seconds=time.time() - t0, r_maxrss_kb_parent=ru.ru_maxrss, staar_converged_log=conv,
                        r_maxrss_note="ru_maxrss of the Rscript parent; the 8 mclapply workers are not summed")
            if rc != 0:
                raise RuntimeError(f"Rscript exited {rc}")
            if conv is not True:
                st = stamp()
                for f in (gf, mf):
                    if f.exists():
                        os.replace(f, f.with_name(f.name + f".rejected_{st}"))
                raise RuntimeError(f"the R log does not contain 'STAAR fit_null_glm: converged TRUE' (found: {conv}); outputs, if any, moved aside with suffix .rejected_{st}")
            ok, rep = validate_outputs(stem, ref_meta, ref_genes, require_converged=True, expect_job_id=my_job)
            meta["validation"] = rep
            if not ok:
                st = stamp()
                for f in (gf, mf):
                    if f.exists():
                        os.replace(f, f.with_name(f.name + f".rejected_{st}"))
                raise RuntimeError(f"outputs failed validation: {rep['problems']} (moved aside with suffix .rejected_{st})")
        meta["r_stage"] = r_stage
        if meta["validation"].get("warnings"):
            print("VALIDATION WARNINGS:", meta["validation"]["warnings"], flush=True)
        if a.compare_existing:
            cmp = compare_with_frozen(stem, a.rep)
            write_json(status_dir / f"{key}_compare.json", cmp)
            meta["compare_with_frozen"] = dict(file=str(status_dir / f"{key}_compare.json"), max_abs_delta_log10p=cmp["max_abs_delta_log10p"], **{"pass": cmp["pass"]})
            print(f"compare with frozen rep{a.rep:03d}: pass={cmp['pass']} max|delta log10 p|={cmp['max_abs_delta_log10p']:.3g}", flush=True)
            if not cmp["pass"]:
                raise RuntimeError("validation gate FAILED: produced p values differ from the frozen result beyond tolerance")
        meta.update(status="accepted", accepted=True, result_genes=str(gf), result_meta=str(mf), result_genes_sha256=sha(gf), result_meta_sha256=sha(mf),
                    result_genes_bytes=gf.stat().st_size)
        rc_final = 0
    except Terminated as e:
        elapsed = time.time() - meta["start"]
        timed_out = STATE["usr1_at"] is not None or (limit_s is not None and elapsed >= limit_s - TIMEOUT_MARGIN)
        meta.update(status="timeout" if timed_out else "terminated", accepted=False, error=str(e), elapsed_at_signal=elapsed,
                    note=("SLURM time limit reached (USR1 warning received or elapsed within 15 min of the limit); counted as an attempt; "
                          "SWEEP.sh resubmits with TIME=12:00:00; two such attempts -> needs_audit") if timed_out else
                         "terminated by signal before the time limit (scancel, node drain or preemption); counted as an attempt")
        print(f"STATUS {key} {meta['status']} after {elapsed / 3600:.2f} h ({e})", flush=True)
        rc_final = 124 if timed_out else 143
    except Exception as e:
        meta.update(status="failed", accepted=False, error=str(e), traceback=traceback.format_exc())
        print(meta["traceback"], flush=True)
        rc_final = 1
    finally:
        meta["elapsed_seconds"] = time.time() - meta["start"]
        meta["finished"] = time.time()
        meta["finished_at"] = now()
        meta["exit_code"] = rc_final
        write_json(status, meta)
        print("STATUS", key, meta["status"], "exit", rc_final, flush=True)
    return rc_final


if __name__ == "__main__":
    sys.exit(main())
