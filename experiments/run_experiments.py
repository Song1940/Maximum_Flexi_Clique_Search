#!/usr/bin/env python3

import argparse, csv, fcntl, gzip, os, re, shutil, signal, statistics, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIELDS = ["group", "dataset", "tau", "algo", "status", "size", "time", "times", "branches",
          "maxrss_kb", "extra", "timeout", "finished_at"]

REAL = "FB GO FW CA CD RO LU DB AM CS WI LJ SK".split()                     # Table II
LFR = [f"lfr_n{n}" for n in (20000, 40000, 60000, 80000, 100000)]            # EQ2, graphs/lfr/
TAU_EQ1 = [0.5, 0.7, 0.9]
TAU_FIG = [k / 100 for k in range(50, 100, 5)]
TAU_LFR = [k / 100 for k in range(90, 100)]
REPEAT_BELOW = 60.0

# algo -> (binary relative to the build directory, --algo)
CPP = {
    "eba": ("flexi", "eba"), "fpa": ("flexi", "fpa"), "npa": ("flexi", "npa"),
    "eba_no_rule1": ("variants/flexi_no_rule1", "eba"), "eba_no_rule2": ("variants/flexi_no_rule2", "eba"),
    "eba_no_rule3": ("variants/flexi_no_rule3", "eba"), "eba_no_incumbent": ("variants/flexi_no_incumbent", "eba"),
    "eba_no_pivot": ("variants/flexi_no_pivot", "eba"), "eba_no_degeneracy": ("variants/flexi_no_degeneracy", "eba"),
    "binary": ("variants/flexi_binary", "eba"), "flat_symse": ("variants/flexi_flat_symse", "eba"),
    "root_stats": ("variants/flexi_root_stats", "eba"),
    "width": ("variants/flexi_width", "eba"), "width_no_pivot": ("variants/flexi_width_no_pivot", "eba"),
}

# group -> algorithms, datasets, taus (width: from the rows of group branching)
GROUPS = {
    "main": dict(algos=["eba", "fpa", "npa"], datasets=REAL, taus=TAU_EQ1),
    "synthetic": dict(algos=["root_stats", "eba"], datasets=LFR, taus=TAU_LFR),
    "effectiveness": dict(algos=["eba", "fpa"], datasets="FB FW CA CD RO LU".split(), taus=TAU_FIG),
    "ablation": dict(algos=["eba", "eba_no_rule1", "eba_no_rule2", "eba_no_rule3", "eba_no_incumbent"],
                     datasets="FB FW CA CD CS WI".split(), taus=TAU_FIG),
    "branching": dict(algos=["eba", "eba_no_pivot", "eba_no_degeneracy", "binary", "flat_symse", "cpsat"],
                      datasets="FB GO FW CA LU CS SK".split(), taus=TAU_FIG),
    "width": dict(algos=["width", "width_no_pivot"], datasets=None, taus=None),
}


def graph_file(d, graph_dir, lfr_dir):
    """Edge list of dataset d; an LFR graph is decompressed next to its .gz file on first use."""
    if not d.startswith("lfr_"):
        return os.path.join(graph_dir, f"{d}.txt")
    txt = os.path.join(lfr_dir, f"{d}.txt")
    if not os.path.exists(txt) and os.path.exists(txt + ".gz"):
        with gzip.open(txt + ".gz", "rb") as fi, open(txt + ".partial", "wb") as fo:
            shutil.copyfileobj(fi, fo)
        os.replace(txt + ".partial", txt)
    return txt


RUNNING = set()                 # process groups of the runs in progress
RUNNING_LOCK = threading.Lock()


def run_once(cmd, timeout, mem_gb):
    """Run cmd with an address-space limit (prlimit) under /usr/bin/time in its own session.
    Returns (stdout, stderr, maxrss_kb, status), status ok / TO / MO / ERR:<last stderr line>."""
    full = ["prlimit", f"--as={mem_gb << 30}", "/usr/bin/time", "-f", "MAXRSS_KB %M"] + cmd
    with RUNNING_LOCK:
        p = subprocess.Popen(full, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                             start_new_session=True)
        RUNNING.add(p.pid)
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.communicate()
        return "", "", "", "TO"
    finally:
        with RUNNING_LOCK:
            RUNNING.discard(p.pid)
    rss = ""
    for line in err.splitlines():
        if line.startswith("MAXRSS_KB"):
            rss = line.split()[1]
    if p.returncode != 0:
        # A failed allocation usually prints one of these; otherwise a failure near the
        # limit is taken as a memory-limit failure too.
        if re.search(r"bad_alloc|MemoryError|Cannot allocate memory|out of memory", err) \
                or (rss and int(rss) >= 0.8 * (mem_gb << 20)):
            return out, err, rss, "MO"
        lines = [l for l in err.strip().splitlines() if not l.startswith("MAXRSS_KB")]
        return out, err, rss, "ERR:" + (lines[-1][:120] if lines else f"exit {p.returncode}").replace(",", ";")
    return out, err, rss, "ok"


def parse_solver(out):
    """'Key : value' lines of the flexi output."""
    kv = {}
    for line in out.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            kv[k.strip()] = v.strip()
    return kv


def solver_extra(algo, err):
    """The 'extra' field: root statistics or the width histogram; None if the line is missing."""
    if algo == "root_stats":
        m = re.search(r"RootStats: F=(\d+) V1=(\d+) V13=(\d+) UB=(\d+)", err)
        return "F={};V1={};V13={};UB={}".format(*m.groups()) if m else None
    if algo in ("width", "width_no_pivot"):
        m = re.search(r"WidthHist:(.*)", err)
        if not m:
            return None
        pairs = [tuple(map(int, x.split(":"))) for x in m.group(1).split()]   # width: count
        internal = sum(c for _, c in pairs)                    # branching nodes (the root level excluded)
        mean = sum(w * c for w, c in pairs) / internal if internal else 0.0
        return (f"internal={internal};mean_w={mean:.6g};max_w={max((w for w, _ in pairs), default=0)};"
                f"capped={sum(c for w, c in pairs if w >= 255)};hist=" + ",".join(f"{w}:{c}" for w, c in pairs))
    return ""


def run_cpp(algo, path, tau, a):
    binary, flag = CPP[algo]
    cmd = [os.path.join(a.build, binary), "--algo", flag, "--file", path, "--tau", str(tau)]
    times, first, failed = [], None, ""
    while True:
        out, err, rss, st = run_once(cmd, a.timeout_of_run, a.mem_gb)
        if st != "ok":
            if first is None:
                return dict(status=st, maxrss_kb=rss)
            failed = st                            # a repeat failed: keep the runs so far, and say so
            break
        kv = parse_solver(out)
        extra = solver_extra(algo, err)
        # A graph without any flexi-clique gives an empty result with Valid "no".
        if (kv.get("Valid") != "yes" and kv.get("Size") != "0") or "Time(s)" not in kv or extra is None:
            return dict(status="ERR:output", maxrss_kb=rss)
        if first is None:
            first = dict(size=kv["Size"], branches=kv.get("Branches", "") if flag == "eba" else "",
                         maxrss_kb=rss, extra=extra)
        elif kv["Size"] != first["size"]:
            return dict(status=f"ERR:repeat size {kv['Size']} != {first['size']}", maxrss_kb=rss)
        times.append(float(kv["Time(s)"]))
        if times[0] >= REPEAT_BELOW or len(times) == 3 or algo in ("width", "width_no_pivot"):
            break                                  # the traces are not timed: one run
    extra = ";".join(x for x in (first.pop("extra"), f"repeat_failed={failed}" if failed else "") if x)
    return dict(first, status="ok", time=f"{statistics.median(times):.6g}",
                times=";".join(f"{x:.6g}" for x in times), extra=extra)


def run_cpsat(path, tau, a):
    cmd = [sys.executable, os.path.join(ROOT, "baselines", "cpsat", "cpsat_baseline.py"), path, str(tau),
           "--fpa", os.path.join(a.build, "flexi"), "--time", str(a.timeout_of_run), "--workers", "1"]
    # CP-SAT stops itself at the time limit; the margin covers building the model.
    out, err, rss, st = run_once(cmd, a.timeout_of_run + 1800, a.mem_gb)
    if st != "ok":
        return dict(status=st, maxrss_kb=rss)
    line = [l for l in out.splitlines() if l.startswith("RESULT")]
    if not line:
        return dict(status="ERR:output", maxrss_kb=rss)
    kv = dict(x.split("=", 1) for x in line[0].split()[1:])
    proved = kv["status"] in ("OPTIMAL", "PROVED_INCUMBENT", "PROVED_BY_PRUNING")
    return dict(status="ok" if proved else "TO", size=kv["size"] if proved else "",
                time=kv["algo"] if proved else "", times=kv["algo"] if proved else "", maxrss_kb=rss,
                extra=f"cpsat_status={kv['status']};Vprime={kv.get('Vprime', '')}")


STAMP = re.compile(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d$")


def read_rows(path):
    """Complete rows of the results file. A row cut off by an interruption (no final newline,
    missing fields, or an incomplete time stamp) is dropped from the file so that appending can
    resume."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        text = f.read()
    header, _, body = text.partition("\n")
    if header.rstrip("\r") != ",".join(FIELDS):
        sys.exit(f"{path}: unexpected header (not a file of this script?)")
    complete = body.endswith("\n") or body == ""
    rows = list(csv.DictReader(body.splitlines(), fieldnames=FIELDS))
    if not complete and rows:
        rows.pop()                                 # the last line has no newline: cut off
    good = [r for r in rows if None not in r and STAMP.match(r.get("finished_at") or "")]
    if len(good) != len(rows) or not complete:
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); w.writerows(good)
    return good


def build_jobs(a, rows):
    sel = lambda opt: {x.strip() for x in opt.split(",") if x.strip()} if opt else None
    want_d, want_a = sel(a.datasets), sel(a.algos)
    want_t = {float(x) for x in sel(a.taus)} if a.taus else None
    jobs = []
    for g in a.groups.split(","):
        spec = GROUPS[g]
        if g == "width":
            insts = sorted({(r["dataset"], float(r["tau"])) for r in rows
                            if r["group"] == "branching" and r["algo"] == "eba" and r["status"] == "ok"
                            and r["branches"] and int(r["branches"]) > 0})
            if not insts:
                print("group width: no EBA row of group branching with branches > 0 yet", flush=True)
        else:
            insts = [(d, t) for d in spec["datasets"] for t in spec["taus"]]
        for d, t in insts:
            for algo in spec["algos"]:
                if (want_d and d not in want_d) or (want_t and t not in want_t) or (want_a and algo not in want_a):
                    continue
                jobs.append((g, d, t, algo))
    return jobs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--groups", required=True, help="comma-separated groups: " + ",".join(GROUPS))
    ap.add_argument("--datasets", help="restrict to these datasets (comma-separated)")
    ap.add_argument("--taus", help="restrict to these values of tau (comma-separated)")
    ap.add_argument("--algos", help="restrict to these algorithms (comma-separated)")
    ap.add_argument("--timeout", type=int, default=86400, help="time limit per run in seconds (default: one day)")
    ap.add_argument("--mem-gb", type=int, default=400, help="memory limit per run in GB (default: 400)")
    ap.add_argument("--workers", type=int, default=1, help="runs in parallel (default: 1)")
    ap.add_argument("--build", default=os.path.join(ROOT, "build"), help="CMake build directory")
    ap.add_argument("--graphs", default=os.path.join(ROOT, "graphs"), help="real-world edge lists")
    ap.add_argument("--lfr", default=os.path.join(ROOT, "graphs", "lfr"), help="LFR graphs")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "results.csv"))
    ap.add_argument("--dry-run", action="store_true", help="list the runs without starting them")
    a = ap.parse_args()
    unknown = [g for g in a.groups.split(",") if g not in GROUPS]
    if unknown:
        ap.error(f"unknown groups: {unknown}")
    if "width" in a.groups.split(",") and "branching" in a.groups.split(","):
        ap.error("group width takes its instances from the rows of group branching: run it in a separate "
                 "invocation after group branching")
    if not shutil.which("prlimit") or not os.path.exists("/usr/bin/time"):
        sys.exit("needs prlimit (util-linux) and /usr/bin/time")

    out = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    lockf = open(out + ".lock", "w")
    try:
        fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)    # one invocation per results file
    except OSError:
        sys.exit(f"another run_experiments.py is writing {out}")
    rows = read_rows(out)
    done = {(r["group"], r["dataset"], float(r["tau"]), r["algo"]) for r in rows if r["status"] in ("ok", "TO", "MO")}
    jobs = [j for j in build_jobs(a, rows) if j not in done]
    print(f"{len(jobs)} runs to do ({len(done)} done in {out})", flush=True)
    if a.dry_run:
        for j in jobs:
            print(*j)
        return

    for d in sorted({j[1] for j in jobs if j[1].startswith("lfr_")}):
        graph_file(d, a.graphs, a.lfr)             # decompress before the workers start
    new = not os.path.exists(out) or os.path.getsize(out) == 0
    fout = open(out, "a", newline="")
    w = csv.DictWriter(fout, fieldnames=FIELDS)
    if new:
        w.writeheader(); fout.flush()
    lock = threading.Lock()
    a.timeout_of_run = a.timeout

    def work(job):
        g, d, t, algo = job
        path = graph_file(d, a.graphs, a.lfr)
        if not os.path.exists(path):
            res = dict(status="ERR:no graph file")
        elif algo == "cpsat":
            res = run_cpsat(path, t, a)
        else:
            res = run_cpp(algo, path, t, a)
        row = dict(group=g, dataset=d, tau=t, algo=algo, timeout=a.timeout,
                   finished_at=time.strftime("%Y-%m-%d %H:%M:%S"), **res)
        with lock:
            w.writerow(row); fout.flush()
            print(f"[{row['finished_at']}] {g} {d} {t} {algo}: {res['status']} "
                  f"size={res.get('size', '')} time={res.get('time', '')}", flush=True)

    def stop(signum, frame):
        # The runs have their own sessions: kill them, and leave without writing their rows.
        with RUNNING_LOCK:
            for pid in RUNNING:
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        os._exit(128 + signum)
    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)

    with ThreadPoolExecutor(max(1, a.workers)) as ex:
        list(ex.map(work, jobs))
    fout.close()


if __name__ == "__main__":
    main()
