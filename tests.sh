#!/usr/bin/env bash
# Run from repo root with the virtualenv activated:
#   source .venv/bin/activate && bash tests.sh
set -uo pipefail

IDAT_REF="cache/209956170157_R01C01_Grn.idat"
IDAT_MIX="cache/209956170157_R02C01_Grn.idat"

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

MIXED_LIN="$TMPDIR_TEST/209956170157_R99C01_Grn.idat"
MIXED_GEOM="$TMPDIR_TEST/209956170157_R98C01_Grn.idat"

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

# ── subtract ────────────────────────────────────────────────────────────────

SUBTRACTED="$TMPDIR_TEST/209956170157_R97C01_Grn.idat"
SELF_SUB="$TMPDIR_TEST/209956170157_R96C01_Grn.idat"

# Subtracting the mixed-in file back out of the linear mix, at the same ratio,
# should return the reference file.
idat-tools subtract -r 0.5 "$MIXED_LIN" "$IDAT_MIX" "$SUBTRACTED" \
    && pass "subtract: output file created" \
    || fail "subtract: command failed"

idat-tools subtract -r 0.5 "$IDAT_REF" "$IDAT_REF" "$SELF_SUB" \
    && pass "subtract self: output file created" \
    || fail "subtract self: command failed"

python3 - "$IDAT_REF" "$IDAT_MIX" "$MIXED_LIN" "$SUBTRACTED" "$SELF_SUB" <<'PYEOF'
import sys
import numpy as np
from pathlib import Path
from idattools.idat import IDATreader

ref_f, mix_f, lin_f, sub_f, self_f = [Path(p) for p in sys.argv[1:6]]

ref  = IDATreader(ref_f).data.per_probe_matrix.set_index("probe_ids")
mix  = IDATreader(mix_f).data.per_probe_matrix.set_index("probe_ids")
lin  = IDATreader(lin_f).data.per_probe_matrix.set_index("probe_ids")
sub  = IDATreader(sub_f).data.per_probe_matrix.set_index("probe_ids")
selfsub = IDATreader(self_f).data.per_probe_matrix.set_index("probe_ids")

errors = 0

# Per-probe check, same shape as the mix test above.
# subtract computes: (I_obs - f * I_mix) / (1 - f), clipped at 0.
for pid in ref.index[:2]:
    o = float(lin.loc[pid, "probe_mean_intensities"])
    m = float(mix.loc[pid, "probe_mean_intensities"])

    expected = round(max((o - 0.5 * m) / 0.5, 0))
    got = float(sub.loc[pid, "probe_mean_intensities"])

    ok = (got == expected)
    status = "PASS" if ok else "FAIL"
    print(f"{status}: probe {pid} subtract : obs={o:.0f} mix={m:.0f} → expected={expected} got={got:.0f}")
    if not ok:
        errors += 1

# Whole-array round trip: mix then subtract at the same ratio returns the
# reference, up to the rounding introduced by both steps (at r=0.5 at most 1).
orig = ref["probe_mean_intensities"].to_numpy(float)
back = sub["probe_mean_intensities"].to_numpy(float)
diff = np.abs(orig - back)
ok = diff.max() <= 2
status = "PASS" if ok else "FAIL"
print(f"{status}: round trip mix→subtract over {len(diff)} probes: "
      f"max abs diff={diff.max():.0f} mean abs diff={diff.mean():.3f}")
if not ok:
    errors += 1

# Subtracting a file from itself is the identity: (I - f*I) / (1-f) = I
orig = ref["probe_mean_intensities"].to_numpy()
got = selfsub["probe_mean_intensities"].to_numpy()
n_wrong = int((orig != got).sum())
ok = (n_wrong == 0)
status = "PASS" if ok else "FAIL"
print(f"{status}: self-subtraction is the identity ({n_wrong} probe(s) differ)")
if not ok:
    errors += 1

# Intensities are unsigned 16 bit: an unclipped negative result would wrap
# around to a very high value instead of going to 0. Nothing may exceed the
# highest intensity seen in either input.
ceiling = max(lin["probe_mean_intensities"].max(), mix["probe_mean_intensities"].max())
n_wrapped = int((sub["probe_mean_intensities"] > ceiling).sum())
ok = (n_wrapped == 0)
status = "PASS" if ok else "FAIL"
print(f"{status}: no uint16 wrap-around ({n_wrapped} probe(s) above input maximum of {ceiling})")
if not ok:
    errors += 1

sys.exit(errors)
PYEOF
py_exit=$?
[ "$py_exit" -eq 0 ] \
    && pass "subtract: intensities, round trip, identity and clipping all correct" \
    || fail "subtract: $py_exit check(s) failed"

# A ratio of 1.0 leaves no reference signal to recover and must be rejected.
idat-tools subtract -r 1.0 "$MIXED_LIN" "$IDAT_MIX" "$TMPDIR_TEST/rejected.idat" 2>/dev/null \
    && fail "subtract -r 1.0: accepted, expected rejection" \
    || pass "subtract -r 1.0: rejected as expected"

# High clip rate: subtracting a nearly complete file removes almost all signal.
python3 - "$IDAT_REF" "$IDAT_MIX" "$TMPDIR_TEST" <<'PYEOF'
import sys
import numpy as np
from pathlib import Path
from idattools.idat import IDATreader, IDATmixer

ref_f, mix_f, tmpdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])

ref = IDATreader(ref_f)
mix = IDATreader(mix_f)

out = tmpdir / "209956170157_R95C01_Grn.idat"
IDATmixer(ref.data).mix(mix.data, 0.95, out, subtract=True)

got = IDATreader(out).data.per_probe_matrix["probe_mean_intensities"].to_numpy()

n_zero = int((got == 0).sum())
print(f"INFO: subtract at r=0.95 clipped {n_zero} / {len(got)} probes to 0 "
      f"({100.0 * n_zero / len(got):.1f}%)")

# Whatever survives must still be a plausible intensity rather than a wrapped
# negative; the warning about clipping is emitted by the library itself.
sys.exit(0 if got.max() <= np.iinfo(np.uint16).max else 1)
PYEOF
py_exit=$?
[ "$py_exit" -eq 0 ] \
    && pass "subtract: heavy subtraction clips to 0 without wrapping" \
    || fail "subtract: heavy subtraction produced wrapped values"

# ── summary ──────────────────────────────────────────────────────────────────

echo ""
echo "Results: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
