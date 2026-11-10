#!/usr/bin/env python3
"""Download the real-world graphs of the paper and convert them to the input format of `flexi`.

  python3 graphs/prepare_datasets.py CS WI LJ            # the graphs that are not included (about 300 MB)
  python3 graphs/prepare_datasets.py                     # all 12 graphs below (the included ones are skipped)
  python3 graphs/prepare_datasets.py --list              # show the graphs and their sources
  python3 graphs/prepare_datasets.py --check build/flexi # check |V| and |E| of every graph (Table II)
  python3 graphs/prepare_datasets.py --communities AM DB # ground-truth communities only

Each graph is written to graphs/<ABBR>.txt as one edge "u<TAB>v" per line, with the node ids
of the original file, after dropping comment lines and self-loops. The files of Network Repository
list the edges in an arbitrary order and some of them twice; these graphs are written in canonical
form (u < v, sorted, each edge once), the KONECT and SNAP graphs in the order of the original file.
These are exactly the files used in the experiments (the order of the edges decides ties in the
search). Downloads are kept in graphs/_download/ and reused. SK-2005 (LAW) needs the
WebGraph library; see graphs/convert_law.sh.

--check runs `flexi --algo fpa` once per graph and compares the size it reports after loading
("Loaded graph: n=... m=...") with |V| and |E| of Table II.

--communities downloads the SNAP ground-truth communities of Amazon and DBLP (Fig. 1, EQ6, EQ7):
graphs/communities/<ABBR>.all.txt (all communities) and <ABBR>.top5000.txt (the 5,000 top-quality
communities), one community per line as tab-separated node ids of graphs/<ABBR>.txt. They
are the decompressed SNAP files, unchanged (same lines, same order), and are checked against the
MD5 sums of the files used in the paper. With --communities, graphs are downloaded only if listed.
Python 3.8+, standard library only.
"""
import argparse, gzip, hashlib, http.client, os, re, shutil, subprocess, sys, tarfile, \
    urllib.request, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
KONECT = "http://konect.cc/files/download.tsv.{0}.tar.bz2"
NR = "https://nrvis.com/download/data/{0}.zip"
SNAP = "https://snap.stanford.edu/data/"
SNAP_COM = SNAP + "bigdata/communities/"

# abbr, name, source, url, member of the archive (None: the compressed file itself), format, |V|, |E|;
# in the order of Table II. Formats: "edges" (Network Repository .edges / .mtx: comma or whitespace
# separated, "%" comments, a non-numeric header line or the MatrixMarket size line is skipped, written in
# canonical form), "konect" and "snap" (whitespace separated, "%" / "#" comments, original order).
GRAPHS = [
    ("FB", "FB-Company", "Network Repository fb-pages-company", NR.format("soc/fb-pages-company"),
     "fb-pages-company.edges", "edges", 14113, 52126),
    ("GO", "Google+", "Network Repository soc-gplus", NR.format("soc/soc-gplus"), "soc-gplus.edges", "edges",
     23628, 39194),
    ("FW", "FrWikinews", "Network Repository ia-frwikinews-user-edits", NR.format("ia/ia-frwikinews-user-edits"),
     "ia-frwikinews-user-edits.edges", "edges", 25042, 68675),
    ("CA", "CAIDA", "KONECT as-caida20071105", KONECT.format("as-caida20071105"),
     "as-caida20071105/out.as-caida20071105", "konect", 26475, 53381),
    ("CD", "Condensed", "KONECT dimacs10-cond-mat-2005", KONECT.format("dimacs10-cond-mat-2005"),
     "dimacs10-cond-mat-2005/out.dimacs10-cond-mat-2005", "konect", 39577, 175691),
    ("RO", "Deezer-RO", "Network Repository soc-gemsec-RO", NR.format("soc/soc-gemsec-RO"), "soc-gemsec-RO.edges",
     "edges", 41773, 125826),
    ("LU", "Luxembourg", "Network Repository inf-luxembourg_osm", NR.format("dimacs10/inf-luxembourg_osm"),
     "inf-luxembourg_osm.mtx", "edges", 114599, 119666),
    ("DB", "DBLP", "SNAP com-DBLP", SNAP_COM + "com-dblp.ungraph.txt.gz", None, "snap", 317080, 1049866),
    ("AM", "Amazon", "SNAP com-Amazon", SNAP_COM + "com-amazon.ungraph.txt.gz", None, "snap", 334863, 925872),
    ("CS", "CiteSeer", "Network Repository co-papers-citeseer", NR.format("dimacs10/co-papers-citeseer"),
     "co-papers-citeseer.mtx", "edges", 434102, 16036720),
    ("WI", "Web-IT", "Network Repository web-it-2004", NR.format("web/web-it-2004"), "web-it-2004.mtx", "edges",
     509338, 7178413),
    ("LJ", "LiveJournal", "SNAP soc-LiveJournal1", SNAP + "soc-LiveJournal1.txt.gz", None, "snap",
     4846609, 42851237),
]
LAW = {"SK": ("SK-2005", 50636059, 1810063330)}
INCLUDED = "FB GO FW CA CD RO LU DB AM".split()            # shipped in graphs/
# SNAP ground-truth communities: abbr -> [(kind, file under SNAP_COM, MD5 of the decompressed file used in the paper)]
COMMUNITIES = {
    "AM": [("all", "com-amazon.all.dedup.cmty.txt.gz", "70b60ac209b2e2b3044b3edbb6195305"),
           ("top5000", "com-amazon.top5000.cmty.txt.gz", "1b5c75db45c1409120b4b1d2854baf82")],
    "DB": [("all", "com-dblp.all.cmty.txt.gz", "d4df167e13599d130a4c5b6e1ed6577c"),
           ("top5000", "com-dblp.top5000.cmty.txt.gz", "4d7c59987314c117183b82107c2bee3a")],
}


def download(url, path):
    """Download url to path, unless path exists; a failed or incomplete download leaves no file."""
    if os.path.exists(path):
        return
    tmp = path + ".part"
    print(f"  downloading {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (flexi-clique dataset script)"})
    got = 0
    try:
        with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
            size = r.headers.get("Content-Length")
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
                got += len(chunk)
        if size is not None and got != int(size):
            raise IOError(f"incomplete download of {url} ({got} of {size} bytes); please run again")
    except (OSError, http.client.HTTPException):     # also URLError, timeouts, truncated reads
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    os.replace(tmp, path)


def open_edges(path, member):
    """Binary stream of the edge file inside the downloaded archive."""
    if path.endswith(".tar.bz2"):
        tar = tarfile.open(path, "r:bz2")
        return tar.extractfile(tar.getmember(member))
    if path.endswith(".zip"):
        return zipfile.ZipFile(path).open(member)
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    raise ValueError(path)


def node_id(tok, line):
    """Node id token without leading zeros; must be a non-negative 32-bit integer (as for flexi)."""
    if not tok.isdigit():
        raise ValueError(f"unexpected line: {line[:80]!r}")
    tok = tok.lstrip(b"0") or b"0"
    if len(tok) >= 10 and int(tok) > 2147483647:
        raise ValueError(f"node id out of the 32-bit range: {line[:80]!r}")
    return tok


def edge_lines(src, fmt):
    """(u, v) byte pairs of the edge lines of src; comments (#, %) and self-loops are dropped. In the
    "edges" format, commas separate the columns of some files, a header line of names is skipped, and
    the first numeric line of a MatrixMarket file (three numbers: rows, columns, entries) is the size."""
    loops = 0
    first = True
    for line in src:
        if fmt == "edges":
            line = line.replace(b",", b" ")
        tok = line.split()
        if not tok or tok[0][:1] in (b"#", b"%"):
            continue
        if fmt == "edges" and first:
            first = False
            if not tok[0].isdigit():                    # header line, e.g. "node_1 node_2"
                continue
            if len(tok) == 3 and line.strip().endswith(tok[2]) and tok[0] == tok[1]:
                continue                                # MatrixMarket size line "n n m"
        if len(tok) < 2:
            raise ValueError(f"unexpected line: {line[:80]!r}")
        u, v = node_id(tok[0], line), node_id(tok[1], line)
        if u == v:
            loops += 1
            continue
        yield u, v
    yield None, loops


def convert(src, dst, header, fmt):
    """Write the edges of src to dst (see edge_lines); canonical form for the "edges" format."""
    tmp = dst + ".part"
    it = edge_lines(src, fmt)
    if fmt == "edges":
        edges = set()
        for u, v in it:
            if u is None:
                loops = v
                break
            a, b = int(u), int(v)
            edges.add((a, b) if a < b else (b, a))
        pairs = ((str(a).encode(), str(b).encode()) for a, b in sorted(edges))
    else:
        pairs = it
        loops = None
    kept = 0
    with src, open(tmp, "wb") as out:
        out.write(header.encode())
        for u, v in pairs:
            if u is None:
                loops = v
                break
            out.write(u + b"\t" + v + b"\n")
            kept += 1
    if kept == 0:
        os.remove(tmp)
        raise ValueError(f"no edges found in the download for {dst}; delete the file in _download and run again")
    os.replace(tmp, dst)
    return kept, loops


def check(binary, path, n, m):
    res = subprocess.run([binary, "--algo", "fpa", "--file", path, "--tau", "0.5"],
                         capture_output=True, text=True)
    if res.returncode != 0:
        return f"FAILED ({binary} exited with status {res.returncode}: {res.stderr.strip()[-200:]})"
    got = re.search(r"Loaded graph: n=(\d+) m=(\d+)", res.stdout)
    if not got:
        return f"FAILED (no 'Loaded graph' line from {binary})"
    gn, gm = int(got.group(1)), int(got.group(2))
    return "OK" if (gn, gm) == (n, m) else f"MISMATCH (n={gn} m={gm}, expected n={n} m={m})"


def prepare_communities(abbr, out):
    """Decompress the SNAP community files of graph abbr into out; return the number of files that
    differ from those used in the paper."""
    dl = os.path.join(out, "_download")
    os.makedirs(dl, exist_ok=True)
    differ = 0
    for kind, fname, md5 in COMMUNITIES[abbr]:
        dst = os.path.join(out, f"{abbr}.{kind}.txt")
        if not os.path.exists(dst):
            arch = os.path.join(dl, fname)
            download(SNAP_COM + fname, arch)
            tmp = dst + ".part"
            with gzip.open(arch, "rb") as src, open(tmp, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
            os.replace(tmp, dst)
        h, lines = hashlib.md5(), 0
        with open(dst, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
                lines += chunk.count(b"\n")
        same = h.hexdigest() == md5
        differ += not same
        print(f"  {dst}: {lines} communities, " + ("identical to the file used in the paper" if same else
              f"DIFFERS from the file used in the paper (MD5 {h.hexdigest()}, expected {md5})"), flush=True)
    return differ


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("graphs", nargs="*", help="abbreviations (default: all)")
    ap.add_argument("--out", default=HERE, help="output directory")
    ap.add_argument("--list", action="store_true", help="list the graphs and exit")
    ap.add_argument("--check", metavar="FLEXI", help="path of the flexi binary; check |V| and |E|")
    ap.add_argument("--communities", nargs="*", metavar="ABBR",
                    help="download the ground-truth communities of these graphs (AM, DB; default both)")
    ap.add_argument("--comm-out", default=os.path.join(HERE, "communities"),
                    help="output directory of the communities")
    args = ap.parse_args()

    if args.list:
        for g in GRAPHS:
            inc = "included" if g[0] in INCLUDED else "download"
            print(f"{g[0]}  {g[1]:<12} |V|={g[6]:>11,} |E|={g[7]:>13,}  {inc:<9} {g[2]}  {g[3]}")
        for a, (name, n, m) in LAW.items():
            print(f"{a}  {name:<12} |V|={n:>11,} |E|={m:>13,}  download  LAW (graphs/convert_law.sh)")
        for a, files in COMMUNITIES.items():
            print(f"{a}  communities: " + ", ".join(f"{a}.{kind}.txt <- {SNAP_COM}{f}" for kind, f, _ in files))
        return
    by_abbr = {g[0]: g for g in GRAPHS}
    comms = None if args.communities is None else [a.upper() for a in args.communities] or list(COMMUNITIES)
    for a in comms or []:
        if a not in COMMUNITIES:
            sys.exit(f"no ground-truth communities for {a}; available: {' '.join(COMMUNITIES)}")
    wanted = [a.upper() for a in args.graphs] or ([] if comms is not None else list(by_abbr))
    for a in wanted:
        if a in LAW:
            sys.exit(f"{a} is a LAW graph; use graphs/convert_law.sh")
        if a not in by_abbr:
            sys.exit(f"unknown graph {a}; see --list")
    if args.check and not os.access(args.check, os.X_OK):
        sys.exit(f"--check: {args.check} is not an executable file")

    differ = failed = 0
    for a in comms or []:
        print(f"[{a}] ground-truth communities", flush=True)
        differ += prepare_communities(a, args.comm_out)
    if wanted:
        os.makedirs(os.path.join(args.out, "_download"), exist_ok=True)
    for a in wanted:
        abbr, name, source, url, member, fmt, n, m = by_abbr[a]
        dst = os.path.join(args.out, abbr + ".txt")
        print(f"[{abbr}] {name}", flush=True)
        if not os.path.exists(dst):
            arch = os.path.join(args.out, "_download", url.rsplit("/", 1)[1])
            download(url, arch)
            header = f"# {name} ({abbr}) from {source}: {url}\n# Table II: |V|={n} |E|={m}\n"
            kept, loops = convert(open_edges(arch, member), dst, header, fmt)
            print(f"  wrote {dst} ({kept} edge lines, {loops} self-loops dropped)")
        else:
            print(f"  {dst} exists")
        if args.check:
            res = check(args.check, dst, n, m)
            failed += res != "OK"
            print(f"  check: {res}", flush=True)
    if failed or differ:
        sys.exit("; ".join(([f"{failed} graph(s) failed the check"] if failed else [])
                           + ([f"{differ} community file(s) differ from those used in the paper"] if differ else [])))


if __name__ == "__main__":
    main()
