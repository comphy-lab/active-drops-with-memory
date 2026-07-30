#!/usr/bin/env python3
import subprocess, sys, shutil

# --------- Settings ----------
EXEC_NAME = "./dropMove"      
SRC       = "dropMove.c"     
AUTO_COMPILE = True

PE_START  = 1.0               # starting guess
STEP0     = 0.5               # initial step size 
TOL       = 0.005             # stop when step < TOL
PE_MIN    = 0.001             # safety floor
PE_MAX    = 100.0             # safety ceiling
MAX_ITERS = 200               # guardrail
# -------------------------------------------

def compile_program():
    cmd = ["qcc", "-O2", "-disable-dimensions", SRC, "-o", EXEC_NAME, "-lm"]
    print("Compiling:", " ".join(cmd))
    subprocess.check_call(cmd)

def run_sim(pe: float) -> bool:
    """
    Runs ./dropMove Pe=<value>. Returns True if MOVED, False if NOT_MOVED.
    """
    out = subprocess.check_output(
        [EXEC_NAME, f"Pe={pe}", "movement_threshold=1"],
                                  universal_newlines=True,
                                  stderr=subprocess.STDOUT)
    moved = None
    # Prefer the last status if multiple lines printed for any reason
    for line in (ln.strip() for ln in out.splitlines() if ln.strip()):
        if line.endswith("STATUS MOVED"):
            moved = True
        elif line.endswith("STATUS NOT_MOVED"):
            moved = False
    if moved is None:
        tail = "\n".join(out.splitlines()[-20:])
        raise RuntimeError(f"Did not find STATUS line in output.\n--- tail ---\n{tail}")
    print(f"Pe={pe:.6g} -> {'MOVED' if moved else 'not moved'}")
    return moved

def clamp(x, lo, hi): return max(lo, min(hi, x))

def find_critical_pe(pe0=PE_START, step0=STEP0):
    if AUTO_COMPILE:
        if not shutil.which("qcc"):
            print("Error: qcc not found; install Basilisk or set AUTO_COMPILE=False.", file=sys.stderr)
            sys.exit(1)
        compile_program()

    pe    = pe0
    step  = step0
    moved = run_sim(pe)

    # Track smallest Pe that still moves (critical edge)
    pe_crit = pe if moved else None

    for _ in range(MAX_ITERS):
        # Move down if moved; up if not moved
        pe_next = clamp(pe - step if moved else pe + step, PE_MIN, PE_MAX)
        moved_next = run_sim(pe_next)

        # Flip in state? halve step and record critical edge
        if moved_next != moved:
            step *= 0.5
            pe_crit = pe_next if moved_next else pe  # first Pe that moves
            print(f"  flip -> step={step:.6g}, provisional Pe_c={pe_crit:.6g}")
            if step < TOL:
                break

        pe, moved = pe_next, moved_next

    if pe_crit is None:
        # never found a moving regime within bounds
        raise RuntimeError("No moving regime found within [PE_MIN, PE_MAX]. Consider widening bounds or increasing STEP0.")

    pe_c_2dp = round(pe_crit, 2)
    print(f"\nCritical Pe (first moving, 2 dp): {pe_c_2dp:.2f}")
    return pe_c_2dp

if __name__ == "__main__":
    # Optional CLI overrides: python find_critical_pe_gd.py [PE_START] [STEP0]
    if len(sys.argv) >= 2:
        PE_START = float(sys.argv[1])
    if len(sys.argv) >= 3:
        STEP0 = float(sys.argv[2])
    find_critical_pe(PE_START, STEP0)
