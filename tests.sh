#!/usr/bin/env bash
# Run from repo root with the virtualenv activated:
#   source .venv/bin/activate && bash tests.sh

set -uo pipefail

IDAT_REF="cache/207513420108_R01C01_Grn.idat"
IDAT_MIX="cache/207513420108_R02C01_Grn.idat"
TMPDIR_TEST="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_TEST"' EXIT

PASS=0
FAIL=0

pass() { echo "PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "FAIL: $1" >&2; FAIL=$((FAIL + 1)); }


# ── view ────────────────────────────────────────────────────────────────────

header_count=$(idat-tools view "$IDAT_REF" | grep -c "^#" || true)
[ "$header_count" -gt 0 ] \
    && pass "view: metadata header present ($header_count lines)" \
    || fail "view: metadata header missing"

no_header_count=$(idat-tools view --no-header "$IDAT_REF" | grep -c "^#" || true)
[ "$no_header_count" -eq 0 ] \
    && pass "view --no-header: no header lines" \
    || fail "view --no-header: header lines still present ($no_header_count)"

n5=$(idat-tools view -n 5   "$IDAT_REF" | grep -cE "^[0-9]" || true)
n15=$(idat-tools view -n 15 "$IDAT_REF" | grep -cE "^[0-9]" || true)
[ "$n5" -lt "$n15" ] \
    && pass "view -n: -n 5 shows fewer rows ($n5) than -n 15 ($n15)" \
    || fail "view -n: expected fewer rows with -n 5 ($n5) than -n 15 ($n15)"


# ── mix ─────────────────────────────────────────────────────────────────────

MIXED_LIN="$TMPDIR_TEST/207513420108_R99C01_Grn.idat"
MIXED_GEOM="$TMPDIR_TEST/207513420108_R98C01_Grn.idat"

idat-tools mix -r 0.5 "$IDAT_REF" "$IDAT_MIX" "$MIXED_LIN" \
    && pass "mix linear: output file created" \
    || fail "mix linear: command failed"

idat-tools mix -r 0.5 --geometric-mean "$IDAT_REF" "$IDAT_MIX" "$MIXED_GEOM" \
    && pass "mix geometric-mean: output file created" \
    || fail "mix geometric-mean: command failed"

python3 - "$IDAT_REF" "$IDAT_MIX" "$MIXED_LIN" "$MIXED_GEOM" <<'PYEOF'
import sys
import numpy as np
from pathlib import Path
from idattools.idat import IDATreader

ref_f, mix_f, lin_f, gm_f = [Path(p) for p in sys.argv[1:5]]

ref = IDATreader(ref_f).data.per_probe_matrix.set_index("probe_ids")
mix = IDATreader(mix_f).data.per_probe_matrix.set_index("probe_ids")
lin = IDATreader(lin_f).data.per_probe_matrix.set_index("probe_ids")
gm  = IDATreader(gm_f).data.per_probe_matrix.set_index("probe_ids")

errors = 0
for pid in ref.index[:2]:
    r = float(ref.loc[pid, "probe_mean_intensities"])
    m = float(mix.loc[pid, "probe_mean_intensities"])

    exp_lin  = round(0.5 * r + 0.5 * m)
    exp_geom = round(np.exp(0.5 * np.log(max(r, 1)) + 0.5 * np.log(max(m, 1))))

    got_lin  = float(lin.loc[pid, "probe_mean_intensities"])
    got_geom = float(gm.loc[pid,  "probe_mean_intensities"])

    for label, expected, got in [
        ("linear   ", exp_lin,  got_lin),
        ("geom-mean", exp_geom, got_geom),
    ]:
        ok = (got == expected)
        status = "PASS" if ok else "FAIL"
        print(f"{status}: probe {pid} {label}: ref={r:.0f} mix={m:.0f} → expected={expected} got={got:.0f}")
        if not ok:
            errors += 1

sys.exit(errors)
PYEOF
py_exit=$?
[ "$py_exit" -eq 0 ] \
    && pass "mix: 2-probe intensity values correct (linear and geometric-mean)" \
    || fail "mix: 2-probe intensity mismatch ($py_exit probe(s) wrong)"


# ── summary ──────────────────────────────────────────────────────────────────

echo ""
echo "Results: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
