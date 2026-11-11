#!/usr/bin/env python3
"""
build_pipeline.py is now a thin front door to build_database.py (the one engine).

Both names work and share the same cache (build_cache/), state, labels and output (build_output/).
Old build_pipeline.py options are translated:

    --genes A B C        ->  --genes A,B,C
    --limit N            ->  --limit N
    --universe all       ->  --include-no-chembl
    --only summarize     ->  --aggregate-only
    --rules-only         ->  --rules-only
    --desc-limit N       ->  --pubchem-desc-limit N
    --refresh-days N     ->  --refresh-older-than N
    --apply-to-db        ->  --merge
    --max-hours H        ->  --max-minutes H*60
    --status             ->  --status
    --work, --verbose    ->  ignored

Any option that build_database.py understands (for example --light-minutes, --heavy-minutes, --light-only)
can be given directly. Files from the earlier pipeline_data/ folder are not read; the first run of this engine
simply fills build_cache/ (downloads that build_database.py already made are reused).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def translate(argv):
    out, i = [], 0
    while i < len(argv):
        a = argv[i]
        if a == "--genes":
            names = []
            i += 1
            while i < len(argv) and not argv[i].startswith("--"):
                names.append(argv[i])
                i += 1
            out += ["--genes", ",".join(names)]
            continue
        if a == "--universe":
            if i + 1 < len(argv) and argv[i + 1] == "all":
                out.append("--include-no-chembl")
            i += 2
            continue
        if a == "--only":
            names = []
            i += 1
            while i < len(argv) and not argv[i].startswith("--"):
                names.append(argv[i])
                i += 1
            if names == ["summarize"]:
                out.append("--aggregate-only")
            continue
        if a == "--desc-limit":
            out += ["--pubchem-desc-limit", argv[i + 1]]
            i += 2
            continue
        if a == "--refresh-days":
            if argv[i + 1] not in ("0", ""):
                out += ["--refresh-older-than", argv[i + 1]]
            i += 2
            continue
        if a == "--apply-to-db":
            out.append("--merge")
        elif a == "--max-hours":
            if float(argv[i + 1]) > 0:
                out += ["--max-minutes", str(float(argv[i + 1]) * 60)]
            i += 2
            continue
        elif a == "--work":
            i += 2
            continue
        elif a == "--verbose":
            pass
        else:
            out.append(a)
        i += 1
    return out


if __name__ == "__main__":
    import build_database
    build_database.main(translate(sys.argv[1:]))
