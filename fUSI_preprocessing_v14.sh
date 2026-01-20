#!/bin/bash

###############################################################################
# fUSI PREPROCESSING PIPELINE — v11 (FULL RESTRUCTURE)
# Steps 0–4 — MATLAB→NIfTI, LocalCopy, Mean, Mask
###############################################################################

# Allow resume
if [[ "$1" == "--continue-from-step" ]]; then
    CONTINUE_FROM="$2"
fi

# Strict mode
set -Eeuo pipefail
trap 'echo -e "\n[ERROR] Failed at line $LINENO" >&2' ERR
export LC_NUMERIC=C
export PYTHONWARNINGS="ignore"
shopt -s extglob

###############################################################################
# COLORS
###############################################################################
RED="\033[38;5;160m"
GREEN="\033[92m"
BLUE="\033[38;5;33m"
PURPLE="\033[38;5;135m"
CYAN="\033[96m"
YELLOW="\033[38;5;220m"
RESET="\033[0m"
BOLD="\033[1m"

###############################################################################
# HELPERS
###############################################################################
center() {
    local msg="$1"
    printf "\n${BOLD}${PURPLE}%*s${RESET}\n" $(( (${#msg}+80)/2 )) "$msg"
}

timestamp(){ date +"%Y%m%d_%H%M%S"; }
msg_info(){ echo -e "${BLUE}[INFO]${RESET} $1"; }
msg_ok(){   echo -e "${GREEN}[OK]${RESET}  $1"; }
msg_skip(){ echo -e "${CYAN}[SKIP]${RESET} $1"; }
msg_warn(){ echo -e "${YELLOW}[WARN]${RESET} $1"; }
msg_fail(){ echo -e "${RED}[FAIL]${RESET} $1"; }

###############################################################################
# GLOBAL LOCAL WORK DIRECTORY (MUST EXIST BEFORE STEP 1)
###############################################################################
WORK="/tmp/fusi_run_$(timestamp)"
mkdir -p "$WORK" || {
    msg_fail "Cannot create local work directory: $WORK"
    exit 1
}

msg_ok "Using local work directory → $WORK"

###############################################################################
# VALIDATION
###############################################################################
if ! command -v python3 >/dev/null; then msg_fail "Python3 missing"; exit 1; fi
if ! command -v fsleyes >/dev/null; then msg_warn "FSLeyes not installed"; fi

###############################################################################
# STEP 0 — SELECT RAW MATLAB FILE
###############################################################################
center "STEP 0 — SELECT RAW MATLAB FILE"
echo -e "${BOLD}Why?${RESET} Load original fUSI MATLAB file (contains I, metadata)."

# -------------------------------------------------------------------------
# Default RAWDATA start directory (GVFS / SMB)
# -------------------------------------------------------------------------
RAWDATA_START="/run/user/3258/gvfs/smb-share:server=wks3,share=pr_ohlendorf/fUS/Project_PACAP_AVATAR_SC/RawData"

# -------------------------------------------------------------------------
# User selects RAW .mat file (dialog OPENS in RawData)
# -------------------------------------------------------------------------
RAW_MAT=$(zenity --file-selection \
    --title="Select fUSI .mat file" \
    --file-filter="MAT files (*.mat) | *.mat" \
    --filename="${RAWDATA_START}/")

[[ -z "$RAW_MAT" ]] && msg_fail "No input file selected." && exit 1

msg_ok "Input MAT file → $RAW_MAT"

# -------------------------------------------------------------------------
# Resolve RawData → AnalysedData mirrored structure
# -------------------------------------------------------------------------
RAW_ROOT=$(dirname "$RAW_MAT")

if [[ "$RAW_ROOT" != *"/RawData/"* ]]; then
    msg_fail "RAW file must be inside a RawData directory"
    exit 1
fi

REL_PATH="${RAW_ROOT#*/RawData/}"
ANALYSED_ROOT="${RAW_ROOT%%/RawData/*}/AnalysedData"
SESSION_ANALYSED="${ANALYSED_ROOT}/${REL_PATH}"

msg_ok "Final analysed folder → $SESSION_ANALYSED"

# -------------------------------------------------------------------------
# Create analysed directory (GVFS-safe)
# -------------------------------------------------------------------------
if mkdir -p "$SESSION_ANALYSED" 2>/dev/null; then
    msg_ok "Analysed directory ready"
else
    msg_warn "Analysed folder is on GVFS / SMB"
    msg_warn "All processing will run in /tmp and results will be copied back"
fi

# -------------------------------------------------------------------------
# Canonical filenames
# -------------------------------------------------------------------------
RAW_FINAL="${SESSION_ANALYSED}/FUS_raw_orient.nii.gz"
PASSTHROUGH_FINAL="${SESSION_ANALYSED}/_passthrough_vars.mat"



###############################################################################
# STEP 1 — MATLAB → NIfTI
###############################################################################
center "STEP 1 — MATLAB → NIfTI"
echo -e "${BOLD}Why?${RESET} Convert MATLAB fUS array into correctly oriented NIfTI."

if [[ -f "$RAW_FINAL" && -f "$PASSTHROUGH_FINAL" ]]; then
    msg_skip "Oriented NIfTI + passthrough vars already exist"
else
    msg_info "Creating oriented NIfTI in /tmp…"

    # ---------------------------------------------------------------------
    # Create WORK directory (ALWAYS defined here)
    # ---------------------------------------------------------------------
    WORK="/tmp/fusi_run_$(timestamp)"
    mkdir -p "$WORK"

    TMP_NII="$WORK/fusi_raw_orient.nii.gz"
    TMP_PASSTHROUGH="$WORK/_passthrough_vars.mat"

python3 <<PY
import numpy as np
import nibabel as nib
import scipy.io as sio
import os

# ------------------------------------------------------------
# Load MATLAB file (FULL schema)
# ------------------------------------------------------------
mat = sio.loadmat(
    "$RAW_MAT",
    squeeze_me=True,
    struct_as_record=False
)

if "I" not in mat:
    raise RuntimeError("MAT file missing required variable 'I'")

I = mat["I"].astype(np.float32)
print("Loaded MATLAB I:", I.shape)

# Expect (Y, X, T)
if I.ndim != 3:
    raise RuntimeError(f"Expected I as (Y,X,T), got {I.shape}")

Y, X, T = I.shape

# ------------------------------------------------------------
# Save ALL passthrough variables (everything except I)
# ------------------------------------------------------------
passthrough = {
    k: v for k, v in mat.items()
    if not k.startswith("__") and k != "I"
}

sio.savemat(
    "$TMP_PASSTHROUGH",
    passthrough,
    do_compression=True
)

print("Saved passthrough vars →", "$TMP_PASSTHROUGH")

# ------------------------------------------------------------
# Orientation fix
# MATLAB: (Y,X,T)
# NIfTI:  (X,1,Y,T)
# + rot90(-1) ensures correct neurological orientation
# ------------------------------------------------------------
I_nifti = np.zeros((X, 1, Y, T), dtype=np.float32)

for t in range(T):
    I_nifti[:, 0, :, t] = np.rot90(I[:, :, t], -1)

# ------------------------------------------------------------
# Affine + voxel sizes
# ------------------------------------------------------------
dx = dy = 0.045
dz = 1.0

aff = np.eye(4)
aff[0,0] = dx
aff[1,1] = dz
aff[2,2] = dy

nii = nib.Nifti1Image(I_nifti, aff)
nii.header.set_zooms((dx, dz, dy, 1.0))

nib.save(nii, "$TMP_NII")
print("Saved oriented NIfTI →", "$TMP_NII")
PY

    # ---------------------------------------------------------------------
    # Copy results BACK to analysed folder (NO mkdir here!)
    # ---------------------------------------------------------------------
    if [[ -d "$SESSION_ANALYSED" ]]; then
        cp "$TMP_NII" "$RAW_FINAL"
        cp "$TMP_PASSTHROUGH" "$PASSTHROUGH_FINAL"
        msg_ok "Copied outputs → $SESSION_ANALYSED"
    else
        msg_warn "Analysed directory unavailable (GVFS)"
        msg_warn "Outputs remain in $WORK"
    fi
fi



###############################################################################
# STEP 2 — LOCAL COPY (FAIL-SAFE)
###############################################################################
center "STEP 2 — LOCAL COPY"
echo -e "${BOLD}Why?${RESET} Avoid slow SMB performance with safe scratch handling."

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PRIMARY_SCRATCH="/tmp"
FALLBACK_SCRATCH="$HOME/scratch"
MIN_FREE_GB=10   # minimum required free space

# Unique run directory
RUN_ID=$(date +%Y%m%d_%H%M%S)
WORK_PRIMARY="${PRIMARY_SCRATCH}/fusi_run_${RUN_ID}"
WORK_FALLBACK="${FALLBACK_SCRATCH}/fusi_run_${RUN_ID}"

# ---------------------------------------------------------------------------
# Function: check free space (GB)
# ---------------------------------------------------------------------------
check_space_gb () {
    df -BG "$1" 2>/dev/null | awk 'NR==2 {gsub("G","",$4); print $4}'
}

# ---------------------------------------------------------------------------
# Decide scratch location
# ---------------------------------------------------------------------------
AVAILABLE_TMP=$(check_space_gb "$PRIMARY_SCRATCH")

if [[ -n "$AVAILABLE_TMP" && "$AVAILABLE_TMP" -ge "$MIN_FREE_GB" ]]; then
    WORK="$WORK_PRIMARY"
    msg_info "Using /tmp scratch (${AVAILABLE_TMP}G free)"
else
    mkdir -p "$FALLBACK_SCRATCH"
    AVAILABLE_FALLBACK=$(check_space_gb "$FALLBACK_SCRATCH")

    if [[ -z "$AVAILABLE_FALLBACK" || "$AVAILABLE_FALLBACK" -lt "$MIN_FREE_GB" ]]; then
        msg_error "Insufficient disk space:"
        msg_error "  /tmp: ${AVAILABLE_TMP:-unknown}G free"
        msg_error "  $FALLBACK_SCRATCH: ${AVAILABLE_FALLBACK:-unknown}G free"
        msg_error "Need at least ${MIN_FREE_GB}G"
        exit 1
    fi

    WORK="$WORK_FALLBACK"
    msg_warn "/tmp full → falling back to $FALLBACK_SCRATCH (${AVAILABLE_FALLBACK}G free)"
fi

mkdir -p "$WORK"

# ---------------------------------------------------------------------------
# Safe copy (atomic)
# ---------------------------------------------------------------------------
RAW_LOCAL="$WORK/FUS_raw_orient.nii.gz"
TMP_COPY="$RAW_LOCAL.tmp"

msg_info "Copying raw data to scratch…"
cp "$RAW_FINAL" "$TMP_COPY" || {
    rm -f "$TMP_COPY"
    msg_error "Copy failed (disk full or I/O error)"
    exit 1
}

mv "$TMP_COPY" "$RAW_LOCAL"

msg_ok "Local copy created → $RAW_LOCAL"

###############################################################################
# STEP 3 — MEAN IMAGE (AUTO-SKIP)
###############################################################################
center "STEP 3 — MEAN IMAGE"
echo -e "${BOLD}Why?${RESET} Required for mask drawing + QC + template building."

MEAN_O="${SESSION_ANALYSED}/mean_FUS_raw_orient.nii.gz"

if [[ -f "$MEAN_O" ]]; then
    msg_skip "Mean image exists → $MEAN_O"
else
    msg_info "Computing mean image…"

python3 <<PY
import nibabel as nib, numpy as np
img = nib.load("$RAW_LOCAL")
d   = img.get_fdata()
mn  = d.mean(axis=-1)
nib.save(nib.Nifti1Image(mn.astype(np.float32), img.affine, img.header), "$MEAN_O")
PY

    msg_ok "Saved → $MEAN_O"
fi


###############################################################################
# STEP 4 — MANUAL MASKING ONLY
###############################################################################
center "STEP 4 — MASK CREATION"
echo -e "${BOLD}Why?${RESET} Brain mask removes edges + noise."

MASK_O="${SESSION_ANALYSED}/mask_FUS_raw_orient.nii.gz"

# -------------------------------------------------------------------------
# If user already saved a mask earlier → skip this step
# -------------------------------------------------------------------------
if [[ -f "$MASK_O" ]]; then
    msg_ok "Mask already exists → $MASK_O"
    msg_ok "Skipping STEP 4."
else

    echo
    echo -e "${BLUE}[INFO]${RESET} Manual masking required."
    echo -e "The mean image will now open in FSLeyes."
    echo -e "Draw brain mask and SAVE AS:"
    echo -e "${RED}$MASK_O${RESET}"
    echo

    # Open MEAN image only (no premask, no auto logic)
    nohup fsleyes "$MEAN_O" >/dev/null 2>&1 &

    # Wait until mask is saved
    while [[ ! -f "$MASK_O" ]]; do
        sleep 2
    done

    msg_ok "Mask saved → $MASK_O"
fi



###############################################################################
# PART 1 COMPLETE
###############################################################################
msg_ok "Part 1/4 complete — ready for Step 5."
###############################################################################
# STEP 5 — CLEAN TIMESERIES (+ FULL MATLAB RESTORE)
###############################################################################
center "STEP 5 — CLEAN TIMESERIES"
echo -e "${BOLD}Why?${RESET} Apply mask and archive clean data for MATLAB."

# -------------------------------------------------------------------------
# FINAL TARGET FILES (SMB — COPY TARGETS)
# -------------------------------------------------------------------------
CLEAN_O="${SESSION_ANALYSED}/cleaned_FUS_raw_orient.nii.gz"
CLEAN_MAT="${SESSION_ANALYSED}/cleaned_FUS_raw_orient.mat"

# -------------------------------------------------------------------------
# TEMP LOCAL FILES (REAL FILESYSTEM)
# -------------------------------------------------------------------------
TMP_CLEAN_NII="${WORK}/cleaned_FUS_raw_orient.nii.gz"
TMP_CLEAN_MAT="${WORK}/cleaned_FUS_raw_orient.mat"

if [[ -f "$CLEAN_O" && -f "$CLEAN_MAT" ]]; then
    msg_skip "Cleaned NIfTI + MAT already exist → skipping Step 5."
else
    msg_info "Applying mask to dataset (local processing)…"

python3 <<PY
import nibabel as nib
import numpy as np
import scipy.io as sio

# ------------------------------------------------------------
# Load NIfTI + mask
# ------------------------------------------------------------
raw  = nib.load("$RAW_LOCAL")
mask = nib.load("$MASK_O").get_fdata() > 0
data = raw.get_fdata().astype(np.float32)

print("Raw NIfTI shape:", data.shape)

# ------------------------------------------------------------
# Apply mask (3D or 4D safe)
# ------------------------------------------------------------
if data.ndim == 4:
    clean = data * mask[..., None]
elif data.ndim == 3:
    clean = data * mask
else:
    raise RuntimeError(f"Unexpected NIfTI shape: {data.shape}")

# ------------------------------------------------------------
# Save CLEANED NIfTI (LOCAL ONLY)
# ------------------------------------------------------------
nib.save(
    nib.Nifti1Image(clean, raw.affine, raw.header),
    "$TMP_CLEAN_NII"
)

# ------------------------------------------------------------
# Convert back to MATLAB orientation (Z,X,T)
# ------------------------------------------------------------
if clean.ndim == 4:
    clean3 = np.squeeze(clean, axis=1)
else:
    clean3 = clean

X, Z, T = clean3.shape

I_mat = np.zeros((Z, X, T), dtype=np.float32)
for t in range(T):
    I_mat[:, :, t] = np.rot90(clean3[:, :, t], +1)

# ------------------------------------------------------------
# Restore passthrough variables
# ------------------------------------------------------------
passthrough = sio.loadmat(
    "$WORK/_passthrough_vars.mat",
    squeeze_me=True,
    struct_as_record=False
)


out = {"I": I_mat}
for k, v in passthrough.items():
    if not k.startswith("__"):
        out[k] = v

# Optional metadata note (non-breaking)
if "metadata" in out:
    try:
        out["metadata"].processing_stage = "Step 5 — masked clean timeseries"
    except Exception:
        pass

# ------------------------------------------------------------
# Save CLEANED MAT (LOCAL ONLY)
# ------------------------------------------------------------
sio.savemat(
    "$TMP_CLEAN_MAT",
    out,
    do_compression=True
)

print("[OK] Local cleaned files written:")
print(" -", "$TMP_CLEAN_NII")
print(" -", "$TMP_CLEAN_MAT")
PY

    # -------------------------------------------------------------------------
    # COPY RESULTS BACK TO SMB / GVFS (BASH ONLY)
    # -------------------------------------------------------------------------
    msg_info "Copying cleaned outputs back to analysed folder (SMB)…"

    mkdir -p "$SESSION_ANALYSED" || {
        msg_fail "Cannot create analysed directory on SMB"
        exit 1
    }

    cp -v "$TMP_CLEAN_NII" "$CLEAN_O" || {
        msg_fail "Failed to copy cleaned NIfTI to analysed folder"
        exit 1
    }

    cp -v "$TMP_CLEAN_MAT" "$CLEAN_MAT" || {
        msg_fail "Failed to copy cleaned MAT to analysed folder"
        exit 1
    }

    msg_ok "Step 5 complete — outputs safely stored on SMB."
fi


###############################################################################
# STEP 6 — ADVANCED QC (AUTO-SKIP IF EXISTS)
###############################################################################
center "STEP 6 — ADVANCED QC"
echo -e "${BOLD}Why?${RESET} Frequency, spatial, temporal, motion-like and stability QC (diagnostic only)."

QC_DIR="${SESSION_ANALYSED}/QC"

if [[ -d "$QC_DIR" ]]; then
    msg_skip "QC folder exists → skipping STEP 6."
else

    mkdir -p "$QC_DIR"

    # Frequency
    FREQ_PNG="$QC_DIR/freq_qc.png"
    FREQ_LF_PNG="$QC_DIR/freq_qc_lowfreq.png"

    # Spatial / temporal
    SPATIAL_PNG="$QC_DIR/spatial_qc.png"
    TEMP_PNG="$QC_DIR/temporal_qc.png"

    # tSNR + SNR/CNR
    TSNR_PNG="$QC_DIR/tsnr_qc.png"
    SNR_PNG="$QC_DIR/snr_cnr_qc.png"

    # Motion-like QC
    MOTION_PNG="$QC_DIR/apparent_motion_xyz.png"

    # Stability QC (Urban/Montaldo)
    STAB_INTENSITY_PNG="$QC_DIR/image_stability_intensity_distribution.png"
    STAB_REJECT_PNG="$QC_DIR/image_stability_rejected_images.png"
    STAB_TRACE="$QC_DIR/image_stability_rejection_trace.png"
    STAB_TXT="$QC_DIR/image_stability_recommendation.txt"

    # Logs
    QC_TXT="$QC_DIR/qc_values.txt"
    WARN_TXT="$QC_DIR/qc_warnings.txt"

    ###############################################################################
    # Window selection (visual + note about indices)
    ###############################################################################
    python3 <<PY
import nibabel as nib, numpy as np, matplotlib.pyplot as plt
img = nib.load("$CLEAN_O")
d = img.get_fdata()
T = d.shape[-1]
gs = d.reshape(-1, T).mean(0)

print("\nIMPORTANT WINDOW SELECTION NOTE")
print("--------------------------------")
print(f"• Total volumes (T) = {T}")
print("• Valid indices: 0 … T-1")
print("• DO NOT enter T as an end index (e.g., if T=3750, last valid is 3749)")
print("• Recommended: leave 1–2 volumes margin at the end\n")

plt.figure(figsize=(12,4))
plt.plot(gs, "k")
plt.title("Global Mean — select QC windows")
plt.xlabel("Volume")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
PY

    read -p "Baseline start end (volumes): " B1 B2
    read -p "Injection start end (volumes): " S1 S2
    read -p "Post-injection start end (volumes): " P1 P2

    ###############################################################################
    # Python QC engine (FAIL-SAFE)  — UNCHANGED
    ###############################################################################
    python3 <<PY
import numpy as np, nibabel as nib, matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq
from scipy.stats import norm

img = nib.load("$CLEAN_O")
d = img.get_fdata()
X,Y,Z,T = d.shape
flat = d.reshape(-1, T)
TR = img.header.get_zooms()[3] if len(img.header.get_zooms())>3 else 1.0

B1,B2 = map(int, ["$B1","$B2"])
S1,S2 = map(int, ["$S1","$S2"])
P1,P2 = map(int, ["$P1","$P2"])

warnings = []

###############################################################################
# 1) Window validation
###############################################################################
def check(name,a,b):
    if a < 0 or b >= T or a >= b:
        warnings.append(f"{name} window invalid: {a}-{b} (T={T})")
        return False
    return True

okB = check("Baseline",B1,B2)
okS = check("Injection",S1,S2)
okP = check("Post",P1,P2)

###############################################################################
# 2) Dead/zero frame detection (helps prevent “post=0” poisoning)
###############################################################################
vol_mean = d.mean(axis=(0,1,2))
med = np.median(vol_mean)
dead = vol_mean < 0.01 * (med + 1e-12)
if dead.any():
    warnings.append(f"{dead.sum()} near-zero frames detected (possible truncation/overwrite).")

###############################################################################
# Helper: FFT
###############################################################################
def fft_block(seg):
    ts = seg.mean(0)
    spec = np.abs(rfft(ts))
    f = rfftfreq(len(ts), TR)
    return f, spec

###############################################################################
# 3) Frequency QC — 0–2 Hz
###############################################################################
plt.figure(figsize=(14,8))
for i,(name,sl,ok) in enumerate([
    ("Baseline",flat[:,B1:B2+1],okB),
    ("Injection",flat[:,S1:S2+1],okS),
    ("Post",flat[:,P1:P2+1],okP),
]):
    plt.subplot(3,1,i+1)
    if not ok:
        plt.text(0.5,0.5,"INVALID WINDOW",ha="center",va="center")
        continue
    f,s = fft_block(sl)
    m = f <= 2
    plt.plot(f[m], s[m])
    plt.title(f"{name} FFT (0–2 Hz)")
    plt.xlabel("Hz")
    plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("$FREQ_PNG", dpi=150); plt.close()

###############################################################################
# 4) Frequency QC — 0–0.1 Hz (normalized, y ticks 0.05)
###############################################################################
plt.figure(figsize=(14,8))
for i,(name,sl,ok) in enumerate([
    ("Baseline",flat[:,B1:B2+1],okB),
    ("Injection",flat[:,S1:S2+1],okS),
    ("Post",flat[:,P1:P2+1],okP),
]):
    plt.subplot(3,1,i+1)
    if not ok:
        plt.text(0.5,0.5,"INVALID WINDOW",ha="center",va="center")
        continue
    f,s = fft_block(sl)
    m = f <= 0.1
    s2 = s[m] / (s[m].max() + 1e-6)
    plt.plot(f[m], s2)
    plt.ylim(0,1.05)
    plt.yticks(np.arange(0,1.01,0.05))
    plt.xlim(0,0.1)
    plt.title(f"{name} FFT (0–0.1 Hz, normalized)")
    plt.xlabel("Hz")
    plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("$FREQ_LF_PNG", dpi=150); plt.close()

###############################################################################
# 5) Spatial QC
###############################################################################
pd_mean = d.mean(-1)
cv = d.std(-1) / (pd_mean + 1e-6)

plt.figure(figsize=(14,8))
plt.subplot(3,1,1); plt.hist(pd_mean.ravel(),80); plt.title("Mean Power Doppler (PD)")
plt.subplot(3,1,2); plt.hist(cv.ravel(),80); plt.title("Temporal CV = std/mean")
plt.subplot(3,1,3); plt.hist(np.diff(d,axis=-1).mean(-1).ravel(),80); plt.title("Axial velocity surrogate (Δt mean)")
plt.tight_layout()
plt.savefig("$SPATIAL_PNG", dpi=150); plt.close()

###############################################################################
# 6) Temporal QC (GS, rGS, DVARS)
###############################################################################
gs = flat.mean(0)
d10 = np.percentile(gs, 10)
rGS = 100*(gs - d10) / (d10 + 1e-6)
DVARS = np.sqrt((np.diff(flat,1,1)**2).mean(0))

plt.figure(figsize=(14,10))
plt.subplot(3,1,1); plt.plot(gs); plt.title("Global signal (GS)")
plt.subplot(3,1,2); plt.plot(rGS); plt.title("Relative GS (rGS) vs 10th percentile baseline")
plt.subplot(3,1,3); plt.plot(DVARS); plt.title("DVARS")
plt.tight_layout()
plt.savefig("$TEMP_PNG", dpi=150); plt.close()

###############################################################################
# 7) tSNR (HEATMAP + HISTOGRAM)
###############################################################################
tsnr = d.mean(-1) / (d.std(-1) + 1e-6)
tsnr_med = float(np.median(tsnr))
mid_y = Y//2

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
im = plt.imshow(tsnr[:,mid_y,:].T, origin="lower", aspect="auto")
plt.colorbar(im)
plt.title("tSNR heatmap (X×Z) at mid-Y")

plt.subplot(1,2,2)
plt.hist(tsnr.ravel(),80)
plt.title(f"tSNR histogram (median={tsnr_med:.2f})")
plt.tight_layout()
plt.savefig("$TSNR_PNG", dpi=150); plt.close()

###############################################################################
# 8) SNR / CNR histogram plot (baseline vs injection)
###############################################################################
if okB:
    baseline = d[...,B1:B2+1]
else:
    baseline = d

if okS:
    inj = d[...,S1:S2+1]
else:
    inj = d

snr = baseline.mean(-1)/(baseline.std(-1)+1e-6)
cnr = np.abs(inj.mean(-1)-baseline.mean(-1)) / (np.sqrt(inj.std(-1)**2 + baseline.std(-1)**2) + 1e-6)

snr_med = float(np.median(snr))
cnr_med = float(np.median(cnr))

plt.figure(figsize=(14,6))
plt.subplot(2,1,1); plt.hist(snr.ravel(),80); plt.title(f"SNR (baseline) median={snr_med:.2f}")
plt.subplot(2,1,2); plt.hist(cnr.ravel(),80); plt.title(f"CNR (inj vs base) median={cnr_med:.2f}")
plt.tight_layout()
plt.savefig("$SNR_PNG", dpi=150); plt.close()

###############################################################################
# 9) Apparent motion QC (center-of-mass drift relative to mid-volume)
###############################################################################
coords = np.indices((X,Y,Z))
ref = d[...,T//2]
ref_sum = np.sum(ref) + 1e-12
c0 = [np.sum(coords[i]*ref)/ref_sum for i in range(3)]

dx,dy,dz = [],[],[]
for t in range(T):
    img_t = d[...,t]
    s = np.sum(img_t) + 1e-12
    c1 = [np.sum(coords[i]*img_t)/s for i in range(3)]
    sh = np.array(c1)-np.array(c0)
    dx.append(sh[0]); dy.append(sh[1]); dz.append(sh[2])

plt.figure(figsize=(12,6))
plt.subplot(3,1,1); plt.plot(dx); plt.title("Apparent drift Δx (voxels) — COM vs mid-volume ref")
plt.subplot(3,1,2); plt.plot(dy); plt.title("Apparent drift Δy (voxels) — COM vs mid-volume ref")
plt.subplot(3,1,3); plt.plot(dz); plt.title("Apparent drift Δz (planes) — COM vs mid-volume ref")
plt.xlabel("Volume")
plt.tight_layout()
plt.savefig("$MOTION_PNG", dpi=150); plt.close()

###############################################################################
# 10) IMAGE STABILITY (Urban/Montaldo) — ROBUST, PAPER-COMPATIBLE
###############################################################################
s = d.mean(axis=(0,1))          # (Z, T)

for iz in range(Z):
    plane = s[iz,:]
    valid = plane[plane > 0]
    if valid.size == 0:
        s[iz,:] = 1.0
        continue
    baseline = np.median(valid[valid < np.percentile(valid, 60)])
    if not np.isfinite(baseline) or baseline == 0:
        baseline = np.median(valid)
    s[iz,:] = plane / (baseline + 1e-12)

lower = s[s <= 1]
if lower.size < 10:
    sigma = 1.4826 * np.median(np.abs(s - 1))
else:
    sigma = 1.4826 * np.median(np.abs(lower - 1))

sigma = max(sigma, 0.02)
threshold = 1 + 3*sigma

outliers = s > threshold
rej_percent = outliers.sum() / outliers.size * 100

plt.figure(figsize=(10,6))
ha, hb = np.histogram(s.ravel(), bins=100)
centers = (hb[:-1] + hb[1:]) / 2
plt.bar(centers, ha, width=(hb[1]-hb[0]), alpha=0.85)

g = np.exp(-0.5*((centers-1)/sigma)**2)
g *= (s.size - outliers.sum()) * (hb[1]-hb[0]) / (sigma*np.sqrt(2*np.pi))
plt.plot(centers, g, 'k', linewidth=2)

plt.axvline(threshold, color='r', linewidth=2)
txt = f"Threshold: {threshold:.3f}\nRejection: {rej_percent:.1f}%"
plt.text(float(threshold), max(ha)*0.6, txt)

plt.title("Intensity distribution (normalized per plane)")
plt.xlabel("Normalized intensity")
plt.ylabel("Number of images (planes × time)")
plt.tight_layout()
plt.savefig("$STAB_INTENSITY_PNG", dpi=150)
plt.close()

t_sec = np.arange(T) * TR
plt.figure(figsize=(12,4))
plt.imshow(1-outliers, aspect='auto', cmap='gray',
           extent=[t_sec[0], t_sec[-1], 0, Z-1], origin='lower')
plt.title("Rejected images (white=accepted, black=rejected)")
plt.xlabel("Time (s)")
plt.ylabel("Planes (Z-slices)")
plt.tight_layout()
plt.savefig("$STAB_REJECT_PNG", dpi=150)
plt.close()

rej_time = outliers.mean(axis=0) * 100
plt.figure(figsize=(12,4))
plt.plot(t_sec, rej_time, 'k', lw=1.2)
plt.axhline(10, color='orange', linestyle='--')
plt.axhline(30, color='red', linestyle='--')
plt.text(t_sec[-1]*0.02, 10.5, "10%: usually acceptable / stable", color='orange')
plt.text(t_sec[-1]*0.02, 30.5, "30%: strong instability likely", color='red')
plt.title("Temporal rejection trace (% rejected planes per frame)")
plt.xlabel("Time (s)")
plt.ylabel("% rejected planes")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("$STAB_TRACE", dpi=150)
plt.close()

total_rej = rej_percent
stim_rej = outliers[:,S1:S2+1].sum()/outliers[:,S1:S2+1].size*100 if okS else 0.0

if total_rej < 10:
    verdict = "STABLE acquisition — no action required."
elif total_rej >= 10 and stim_rej < total_rej*0.5:
    verdict = ("MECHANICAL instability suspected.\n"
               "Check headpost, probe fixation, and coupling.")
else:
    verdict = ("STIMULUS-LOCKED artifacts detected.\n"
               "Consider a blind time window around stimulation.\n"
               "Averaging cannot recover this interval.")

with open("$STAB_TXT","w") as f:
    f.write(f"Total rejected frames: {total_rej:.2f}%\n")
    f.write(f"Stimulus-window rejected frames: {stim_rej:.2f}%\n")
    f.write(f"Threshold: {threshold:.4f}\n")
    f.write(f"Sigma: {sigma:.4f}\n\n")
    f.write("QC Recommendation:\n")
    f.write(verdict + "\n")

###############################################################################
# 11) Save QC values + warnings
###############################################################################
with open("$QC_TXT","w") as f:
    f.write(f"Volumes (T): {T}\n")
    f.write(f"TR (s): {TR}\n")
    f.write(f"tSNR median: {tsnr_med:.4f}\n")
    f.write(f"SNR median: {snr_med:.4f}\n")
    f.write(f"CNR median: {cnr_med:.4f}\n")
    f.write(f"Stability sigma: {sigma:.6f}\n")
    f.write(f"Stability threshold: {threshold:.6f}\n")
    f.write(f"Total rejected (%): {total_rej:.2f}\n")

if warnings:
    with open("$WARN_TXT","w") as f:
        for w in warnings:
            f.write("[WARNING] "+w+"\n")

print("STEP 6 QC complete")
PY

    msg_ok "STEP 6 complete → $QC_DIR"

fi



###############################################################################
# STEP 7 — TRIMMING
###############################################################################
center "STEP 7 — TRIMMING"
echo -e "${BOLD}Why?${RESET} Remove unstable initial/final periods."

TRIM_O="$CLEAN_O"

read -p "Skip trimming? (y/n): " SKIP_TRIM
if [[ "$SKIP_TRIM" != "y" ]]; then
    echo -e "${BLUE}Select file for trimming:${RESET}"

    mapfile -t FILES < <(find "$SESSION_ANALYSED" -maxdepth 1 -type f -name "*.nii.gz" | sort)
    for i in "${!FILES[@]}"; do printf "%3d) %s\n" $((i+1)) "$(basename "${FILES[$i]}")"; done

    read -p "Choice: " IDX
    SELFILE="${FILES[$((IDX-1))]}"

    TR=$(python3 - <<EOF
import nibabel as nib
z=nib.load("$SELFILE").header.get_zooms()
print(1.0 if len(z)<4 else z[3])
EOF
)
  

    read -p "Trim START seconds: " TS_SEC
    read -p "Trim END seconds: " TE_SEC

python3 <<PY
import nibabel as nib, numpy as np

img=nib.load("$SELFILE")
d=img.get_fdata()
TR=float("$TR")
s=int(float("$TS_SEC")/TR)
e=int(float("$TE_SEC")/TR)

d2 = d[..., s:d.shape[-1]-e]
out="$SESSION_ANALYSED/trim_$(timestamp).nii.gz"
nib.save(nib.Nifti1Image(d2.astype(np.float32), img.affine,img.header),out)
print(out)
PY

    TRIM_O=$(tail -1 <<< "$(python3 - <<EOF
print("")
EOF
)")
fi

###############################################################################
# STEP 8 — SMOOTHING
###############################################################################
center "STEP 8 — SMOOTHING"
echo -e "${BOLD}Why?${RESET} Temporal + spatial smoothing increases SNR."

SMOOTH_DIR="${SESSION_ANALYSED}/SmoothedData"
mkdir -p "$SMOOTH_DIR"

SMOOTH_O="$TRIM_O"

read -p "Skip smoothing? (y/n): " SKIP_SMOOTH
if [[ "$SKIP_SMOOTH" != "y" ]]; then

    echo -e "${BLUE}Select file for smoothing:${RESET}"
    mapfile -t FILES < <(find "$SESSION_ANALYSED" -maxdepth 1 -type f -name "*.nii.gz" | sort)
    for i in "${!FILES[@]}"; do printf "%3d) %s\n" $((i+1)) "$(basename "${FILES[$i]}")"; done

    read -p "Choice: " IDX
    SELFILE="${FILES[$((IDX-1))]}"

    TR=$(python3 - <<EOF
import nibabel as nib
z=nib.load("$SELFILE").header.get_zooms()
print(1.0 if len(z)<4 else z[3])
EOF
)

    ###############################################
    # TEMPORAL SMOOTHING
    ###############################################
    read -p "Temporal smoothing window (seconds, 0=skip): " TSM
    if [[ "$TSM" != "0" ]]; then

        TMP="${SMOOTH_DIR}/ts_${TSM}s_$(timestamp).nii.gz"

python3 <<PY
import nibabel as nib, numpy as np
from scipy.ndimage import uniform_filter1d

img=nib.load("$SELFILE")
d=img.get_fdata()
win=max(1,int(float("$TSM")/float("$TR")))
out = uniform_filter1d(d, size=win, axis=-1, mode="nearest")

nib.save(nib.Nifti1Image(out.astype(np.float32),img.affine,img.header),"$TMP")
PY

        msg_ok "Temporal smoothing → $TMP"
        SELFILE="$TMP"
    fi

    ###############################################
    # SPATIAL SMOOTHING
    ###############################################
    read -p "Spatial smoothing FWHM (mm, 0=skip): " FWHM
    if [[ "$FWHM" != "0" ]]; then

        SIG=$(python3 - <<EOF
import numpy as np
print(float("$FWHM")/np.sqrt(8*np.log(2)))
EOF
)

        TMP="${SMOOTH_DIR}/ss_${FWHM}mm_$(timestamp).nii.gz"

python3 <<PY
import nibabel as nib, numpy as np
from scipy.ndimage import gaussian_filter

img=nib.load("$SELFILE")
d=img.get_fdata()
sigma=float("$SIG")

out=np.stack([gaussian_filter(d[...,t],sigma)
              for t in range(d.shape[-1])],axis=-1)

nib.save(nib.Nifti1Image(out.astype(np.float32),img.affine,img.header),"$TMP")
PY

        SMOOTH_O="$TMP"
    else
        SMOOTH_O="$SELFILE"
    fi
fi

msg_ok "Smoothing output → $SMOOTH_O"

###############################################################################
# END PART 2
###############################################################################
msg_ok "Part 2/4 complete — ready for scrubbing."

###############################################################################
# STEP 9 — SCRUBBING (rGS / rDVARS strict thresholds)
###############################################################################
center "STEP 9 — SCRUBBING (OPTIONAL)"
echo -e "${BOLD}Why?${RESET} Remove motion-contaminated frames using strict QC thresholds."
echo -e "${YELLOW}NOTE:${RESET} Use this INSTEAD of frame-rate rejection (Step 6), not in addition."

SCRUB_DIR="${SESSION_ANALYSED}/ScrubbedData"
mkdir -p "$SCRUB_DIR"

read -p "Skip scrubbing? (y/n): " SKIP_SCRUB
if [[ "$SKIP_SCRUB" != "y" ]]; then

    ###########################################################################
    # FILE SELECTION
    ###########################################################################
    echo -e "${BLUE}Select file for scrubbing:${RESET}"
    mapfile -t FILES < <(find "$SESSION_ANALYSED" -type f -name "*.nii.gz" | sort)

    for i in "${!FILES[@]}"; do
        printf "%3d) %s\n" $((i+1)) "${FILES[$i]#$SESSION_ANALYSED/}"
    done

    read -p "Choice: " IDX
    SELFILE="${FILES[$((IDX-1))]}"

    SCRUB_O="$SCRUB_DIR/scrubbed_rGS_rDVARS_$(timestamp).nii.gz"

###############################################################################
# PYTHON SCRUBBING + QC
###############################################################################
python3 <<PY
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import os, sys

img = nib.load("$SELFILE")
d = img.get_fdata()
T = d.shape[-1]
flat = d.reshape(-1, T)

###############################################################################
# Metrics BEFORE
###############################################################################
gs = flat.mean(axis=0)
D1 = np.percentile(gs, 10)
rGS = 100 * (gs - D1) / (D1 + 1e-6)

diff = np.diff(flat, axis=1)
DVARS = np.sqrt((diff**2).mean(axis=0))
DVARS_full = np.concatenate(([0], DVARS))

p25, p75 = np.percentile(flat, [25, 75], axis=1)
iqr = p75 - p25
sigma = iqr / 1.349
mu0 = np.sqrt((2 * sigma**2).mean())
rDVARS = DVARS_full / (mu0 + 1e-6)

###############################################################################
# Thresholds
###############################################################################
TH_RGS = 10.0
TH_RDV = 2.0

bad = np.where((rGS > TH_RGS) | (rDVARS > TH_RDV))[0]
good = np.setdiff1d(np.arange(T), bad)

###############################################################################
# Replace bad frames
###############################################################################
d2 = d.copy()
mean_good = d[..., good].mean(axis=-1)
for b in bad:
    d2[..., b] = mean_good

nib.save(nib.Nifti1Image(d2.astype(np.float32), img.affine, img.header), "$SCRUB_O")

###############################################################################
# Metrics AFTER
###############################################################################
flat2 = d2.reshape(-1, T)
gs2 = flat2.mean(axis=0)
diff2 = np.diff(flat2, axis=1)
DVARS2 = np.concatenate(([0], np.sqrt((diff2**2).mean(axis=0))))

###############################################################################
# QC OUTPUTS
###############################################################################
qcdir = "$SESSION_ANALYSED/QC"
os.makedirs(qcdir, exist_ok=True)
tag = os.path.basename("$SCRUB_O").replace(".nii.gz","")

# 1) Metrics + rejected frames
plt.figure(figsize=(14,10))
plt.subplot(3,1,1)
plt.plot(rGS, 'b'); plt.axhline(TH_RGS, c='r', ls='--')
plt.scatter(bad, rGS[bad], c='r', s=20)
plt.title("rGS with scrubbed volumes")

plt.subplot(3,1,2)
plt.plot(rDVARS, 'k'); plt.axhline(TH_RDV, c='r', ls='--')
plt.scatter(bad, rDVARS[bad], c='r', s=20)
plt.title("rDVARS with scrubbed volumes")

plt.subplot(3,1,3)
mask = np.zeros(T); mask[bad] = 1
plt.imshow(mask[None,:], aspect="auto", cmap="Reds")
plt.yticks([])
plt.xlabel("Time (volumes)")
plt.title("Scrub mask (red = replaced)")

plt.tight_layout()
plt.savefig(f"{qcdir}/scrubbing_detection_{tag}.png", dpi=150)
plt.close()

# 2) Before vs After GS
plt.figure(figsize=(12,4))
plt.plot(gs, label="Before", alpha=0.7)
plt.plot(gs2, label="After", alpha=0.7)
plt.legend(); plt.title("Global Signal: Before vs After Scrubbing")
plt.tight_layout()
plt.savefig(f"{qcdir}/scrubbing_GS_before_after_{tag}.png", dpi=150)
plt.close()

# 3) Before vs After DVARS
plt.figure(figsize=(12,4))
plt.plot(DVARS_full, label="Before", alpha=0.7)
plt.plot(DVARS2, label="After", alpha=0.7)
plt.legend(); plt.title("DVARS: Before vs After Scrubbing")
plt.tight_layout()
plt.savefig(f"{qcdir}/scrubbing_DVARS_before_after_{tag}.png", dpi=150)
plt.close()

print("[OK] Scrubbing QC saved in:", qcdir)
PY

    msg_ok "Scrubbing complete → $SCRUB_O"

else
    msg_skip "Scrubbing skipped."
    SCRUB_O="$SMOOTH_O"
fi



###############################################################################
# STEP 10 — FILTERING (Butterworth / SVD / eSVD)
###############################################################################
center "STEP 10 — FILTERING (Butterworth / SVD / eSVD)"
echo -e "${BOLD}Why?${RESET} Remove tissue clutter (SVD / eSVD) and preserve hemodynamic frequencies."

FILT_DIR="${SESSION_ANALYSED}/FilteredData"
mkdir -p "$FILT_DIR"

QC_DIR="${SESSION_ANALYSED}/QC"
mkdir -p "$QC_DIR"

step10_loop=true
while $step10_loop; do

    # --------------------------------------------------------
    # SKIP?
    # --------------------------------------------------------
    read -p "Skip filtering? (y/n): " SKIP_FILT
    if [[ "$SKIP_FILT" == "y" ]]; then
        msg_skip "Filtering skipped."
        # Use previous best dataset as input for Step 11
        FILT_O="${SCRUB_O:-$SMOOTH_O}"
        break
    fi

    # --------------------------------------------------------
    # FILE SELECTION
    # --------------------------------------------------------
    echo "Select file for filtering:"
    mapfile -t FILES < <(find "$SESSION_ANALYSED" -type f -name "*.nii.gz" | sort)

    if [[ "${#FILES[@]}" -eq 0 ]]; then
        msg_fail "No NIfTI files found in $SESSION_ANALYSED"
        break
    fi

    for i in "${!FILES[@]}"; do
        printf "%3d) %s\n" "$((i+1))" "$(basename "${FILES[$i]}")"
    done

    read -p "Choice: " IDX
    SELFILE="${FILES[$((IDX-1))]}"
    SELNAME=$(basename "$SELFILE")
    TIMESTAMP=$(timestamp)

    msg_info "Selected file → $SELFILE"
    echo

    # --------------------------------------------------------
    # TR DETECTION + OVERRIDE
    # --------------------------------------------------------
    TR_META=$(python3 - <<EOF
import nibabel as nib
img = nib.load("$SELFILE")
z = img.header.get_zooms()
print(z[3] if len(z) >= 4 else -1.0)
EOF
)
    echo -e "${BLUE}Detected TR from header: ${TR_META} s${RESET}"
    read -p "Press ENTER to accept, or type TR in seconds: " TR_INPUT
    if [[ -z "$TR_INPUT" ]]; then
        TR="$TR_META"
    else
        TR="$TR_INPUT"
    fi
    msg_ok "Using TR = $TR s"

    # --------------------------------------------------------
    # FILTER TYPE MENU
    # --------------------------------------------------------
    echo -e "${BLUE}Choose filter type:${RESET}"
    echo "1) High-pass (Butterworth)"
    echo "2) Low-pass (Butterworth)"
    echo "3) Band-pass (Butterworth)"
    echo "4) Classical SVD clutter filtering"
    echo "5) SVD + High-pass hybrid"
    echo "6) Enhanced SVD (eSVD for magnitude data)"
    read -p "Choice: " FTYPE

    HP=0; LP=0; BP1=0; BP2=0

    case "$FTYPE" in
        1)
            read -p "High-pass cutoff (Hz): " HP
            ;;
        2)
            read -p "Low-pass cutoff (Hz): " LP
            ;;
        3)
            read -p "Band-pass LOW cutoff (Hz): " BP1
            read -p "Band-pass HIGH cutoff (Hz): " BP2
            ;;
        4|5|6)
            echo -e "${PURPLE}[INFO] SVD-based filtering selected — SVD QC will be generated.${RESET}"
            ;;
        *)
            msg_warn "Invalid filter type; skipping Step 10."
            FILT_O="$SELFILE"
            break
            ;;
    esac

    FILT_O="${FILT_DIR}/filt_${FTYPE}_${TIMESTAMP}.nii.gz"

    # --------------------------------------------------------
    # PYTHON ENGINE (Butterworth / SVD / eSVD)
    # --------------------------------------------------------
python3 <<PY
import numpy as np, nibabel as nib, matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
import os

fname     = "$SELFILE"
outname   = "$FILT_O"
ftype     = "$FTYPE"
tr        = float("$TR")
nyq       = 0.5 / tr
qcdir     = "$QC_DIR"
timestamp = "$TIMESTAMP"
selname   = "$SELNAME"

HP  = float("$HP")
LP  = float("$LP")
BP1 = float("$BP1")
BP2 = float("$BP2")

os.makedirs(qcdir, exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
img = nib.load(fname)
d   = img.get_fdata().astype(np.float32)
aff = img.affine
hdr = img.header
X, Y, Z, T = d.shape
flat = d.reshape(-1, T)

# ----------------------------
# GLOBAL MEAN BEFORE FILTER
# ----------------------------
gs_before = flat.mean(0)
plt.figure(figsize=(12,4))
plt.plot(gs_before, "k")
plt.title(f"Global Mean BEFORE filtering ({selname})")
plt.xlabel("Volume")
plt.grid(alpha=0.3)
plt.tight_layout()
gm_before_png = os.path.join(qcdir, f"GM_before_{ftype}_{selname}_{timestamp}.png")
plt.savefig(gm_before_png, dpi=150)
plt.close()
print("[QC] Saved:", gm_before_png)

# ----------------------------
# FFT BEFORE FILTER
# ----------------------------
freqs = np.fft.rfftfreq(T, tr)
fft_before = np.abs(np.fft.rfft(gs_before))

plt.figure(figsize=(10,5))
plt.plot(freqs, fft_before)
plt.title(f"FFT BEFORE filtering ({selname})")
plt.xlabel("Hz")
plt.grid(alpha=0.3)
plt.tight_layout()
fft_before_png = os.path.join(qcdir, f"FFT_before_{ftype}_{selname}_{timestamp}.png")
plt.savefig(fft_before_png, dpi=150)
plt.close()
print("[QC] Saved:", fft_before_png)

# ----------------------------
# SVD helper: spectrum + auto elbow
# ----------------------------
def svd_with_qc(F):
    print("[INFO] Computing SVD for QC...")
    U, s, Vt = np.linalg.svd(F, full_matrices=False)

    # Linear spectrum
    plt.figure(figsize=(6,4))
    plt.plot(s, "k")
    plt.title("SVD spectrum (linear)")
    plt.xlabel("Component")
    plt.ylabel("Singular value")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    svd_lin_png = os.path.join(qcdir, f"SVD_spectrum_linear_{selname}_{timestamp}.png")
    plt.savefig(svd_lin_png, dpi=150)
    plt.close()
    print("[QC] Saved:", svd_lin_png)

    # Log spectrum
    plt.figure(figsize=(6,4))
    plt.semilogy(s, "k")
    plt.title("SVD spectrum (log)")
    plt.xlabel("Component")
    plt.ylabel("Singular value")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    svd_log_png = os.path.join(qcdir, f"SVD_spectrum_log_{selname}_{timestamp}.png")
    plt.savefig(svd_log_png, dpi=150)
    plt.close()
    print("[QC] Saved:", svd_log_png)

    # Simple elbow detector: max drop
    ds = np.abs(np.diff(s))
    elbow = int(np.argmax(ds)) + 1
    elbow = max(1, elbow)

    print(f"[INFO] Auto tissue cutoff N (elbow) = {elbow}")
    # Save this info to text file
    with open(os.path.join(qcdir, f"SVD_params_{selname}_{timestamp}.txt"), "w") as f:
        f.write(f"Auto tissue cutoff (elbow) N = {elbow}\n")
        f.write(f"First 20 singular values: {s[:20]}\n")

    return U, s, Vt, elbow

# ----------------------------
# Butterworth filter
# ----------------------------
def butter_filter(F):
    if ftype == "1":
        b, a = butter(6, HP/nyq, btype="highpass")
    elif ftype == "2":
        b, a = butter(6, LP/nyq, btype="lowpass")
    elif ftype == "3":
        b, a = butter(6, [BP1/nyq, BP2/nyq], btype="bandpass")
    else:
        return F

    out = np.zeros_like(F)
    for i in range(F.shape[0]):
        x = F[i]
        if np.isfinite(x).all() and x.std() > 1e-8:
            try:
                out[i] = filtfilt(b, a, x)
            except Exception:
                out[i] = x
        else:
            out[i] = x
    return out

# ----------------------------
# eSVD for magnitude-only data
# ----------------------------
def eSVD_filter(F):
    U, s, Vt, elbow = svd_with_qc(F)

    # Remove tissue components (first N_tissue)
    N_tissue = elbow
    s2 = s.copy()
    s2[:N_tissue] = 0.0

    blood_idx = np.where(s2 > 0)[0]
    if blood_idx.size == 0:
        print("[WARN] No blood subspace found; returning original data.")
        return F

    # Partition blood components into K subspaces with ~equal energy
    K = 5
    energies = (s2[blood_idx] ** 2)
    total_E = energies.sum()
    target = total_E / K
    subspaces = []
    current = []
    acc = 0.0

    for idx, e in zip(blood_idx, energies):
        current.append(idx)
        acc += e
        if acc >= target and len(subspaces) < K - 1:
            subspaces.append(np.array(current))
            current = []
            acc = 0.0
    if current:
        subspaces.append(np.array(current))

    U_w = U.copy()

    # For each blood subspace, compute feature map and weight U
    for sub in subspaces:
        if sub.size == 0:
            continue
        # Weighted spatial singular images (approx blood content)
        Asub = (U[:, sub] * s2[sub])  # (Nvox, n_sub)
        proj = np.mean(np.abs(Asub), axis=1)  # (Nvox,)

        proj_img = proj.reshape(X, Y, Z)
        proj_norm = proj_img / (proj_img.max() + 1e-6)

        # Gentle enhancement: gamma < 1 boosts weaker vessel signals
        gamma = 0.7
        w = proj_norm ** gamma  # (X,Y,Z)

        # Apply weights back to U columns
        w_flat = w.reshape(-1, 1)
        U_w[:, sub] = U_w[:, sub] * w_flat

    F_filt = (U_w * s2) @ Vt
    return F_filt

# ----------------------------
# APPLY FILTER TYPE
# ----------------------------
if ftype in ["1", "2", "3"]:
    flat_f = butter_filter(flat)

elif ftype == "4":
    # Classical SVD clutter filter with auto cutoff
    U, s, Vt, N_tissue = svd_with_qc(flat)
    s2 = s.copy()
    s2[:N_tissue] = 0.0
    flat_f = (U * s2) @ Vt

elif ftype == "5":
    # SVD + high-pass hybrid
    U, s, Vt, N_tissue = svd_with_qc(flat)
    s2 = s.copy()
    s2[:N_tissue] = 0.0
    S_svd = (U * s2) @ Vt
    flat_f = butter_filter(S_svd)

elif ftype == "6":
    flat_f = eSVD_filter(flat)

else:
    print("[WARN] Invalid filter type; copying input.")
    flat_f = flat.copy()

# ----------------------------
# SAVE FILTERED DATA
# ----------------------------
d_f = flat_f.reshape(X, Y, Z, T).astype(np.float32)
nib.save(nib.Nifti1Image(d_f, aff, hdr), outname)
print("[OK] Saved filtered dataset →", outname)

# ----------------------------
# GLOBAL MEAN AFTER FILTER
# ----------------------------
gs_after = flat_f.mean(0)
plt.figure(figsize=(12,4))
plt.plot(gs_after, "b")
plt.title(f"Global Mean AFTER filtering ({selname})")
plt.xlabel("Volume")
plt.grid(alpha=0.3)
plt.tight_layout()
gm_after_png = os.path.join(qcdir, f"GM_after_{ftype}_{selname}_{timestamp}.png")
plt.savefig(gm_after_png, dpi=150)
plt.close()
print("[QC] Saved:", gm_after_png)

# ----------------------------
# FFT AFTER FILTER + OVERLAY
# ----------------------------
fft_after = np.abs(np.fft.rfft(gs_after))

plt.figure(figsize=(10,5))
plt.plot(freqs, fft_before, "r", label="Before")
plt.plot(freqs, fft_after,  "b", label="After")
plt.title(f"FFT Before vs After filtering ({selname})")
plt.xlabel("Hz")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
fft_ba_png = os.path.join(qcdir, f"FFT_before_after_{ftype}_{selname}_{timestamp}.png")
plt.savefig(fft_ba_png, dpi=150)
plt.close()
print("[QC] Saved:", fft_ba_png)

PY

    msg_ok "Step 10 finished → $FILT_O"

    # ----------------------------------------------------
    # REPEAT STEP 10?
    # ----------------------------------------------------
    read -p "Repeat Step 10 with another file? (y/n): " REP
    if [[ "$REP" != "y" ]]; then
        step10_loop=false
    fi

done

###############################################################################
# STEP 11 — CONFOUND REGRESSION + OPTIONAL PCA
###############################################################################
center "STEP 11 — CONFOUND REGRESSION + PCA"
echo -e "${BOLD}Why?${RESET} Remove global/physiological/motion-like confounds (GSR/CompCor), \
then optionally apply PCA denoising."

CR_DIR="${SESSION_ANALYSED}/ConfoundRegression"
mkdir -p "$CR_DIR"

read -p "Skip Step 11 (confound regression + PCA)? (y/n): " SKIP_11
if [[ "$SKIP_11" == "y" ]]; then
    msg_skip "Skipping Step 11."
else

###############################################################################
# SELECT INPUT FILE FOR CONFOUND REGRESSION
###############################################################################
echo -e "${BLUE}Select file for confound regression:${RESET}"

mapfile -t FILES < <(find "$SESSION_ANALYSED" -type f -name "*.nii.gz" | sort)
for i in "${!FILES[@]}"; do
    printf "%3d) %s\n" $((i+1)) "$(basename "${FILES[$i]}")"
done

read -p "Choice: " IDX
SELFILE="${FILES[$((IDX-1))]}"
SELNAME=$(basename "$SELFILE")
msg_ok "Selected → $SELNAME"

###############################################################################
# SELECT CONFOUND METHOD
###############################################################################
echo -e "${BLUE}Select confound regression method:${RESET}"
echo "1) None"
echo "2) Global Signal Regression (GSR)"
echo "3) tCompCor (5 comps)"
echo "4) aCompCor (WM+CSF, 5 comps)"
echo "5) Random CompCor (5 comps)"
echo "6) Low-variance CompCor (5 comps)"
echo "7) PCA only (no confound regression)"
read -p "Choice: " CTYPE

timestamp_now=$(date +"%Y%m%d_%H%M%S")
CONF_O="$SELFILE"

###############################################################################
# BUILD & APPLY CONFOUND REGRESSORS
###############################################################################
if [[ "$CTYPE" != "1" && "$CTYPE" != "7" ]]; then
    CONF_O="${CR_DIR}/confreg_${CTYPE}_${timestamp_now}.nii.gz"

    # aCompCor requires WM+CSF mask
    if [[ "$CTYPE" == "4" ]]; then
        echo -e "${PURPLE}[ACTION REQUIRED]${RESET}"
        echo "Open mean image: $MEAN_O"
        echo "Draw WM mask → save as $SESSION_ANALYSED/wm_mask.nii.gz"
        echo "Draw CSF mask → save as $SESSION_ANALYSED/csf_mask.nii.gz"
        echo "Waiting..."
        while [[ ! -f "$SESSION_ANALYSED/wm_mask.nii.gz" || ! -f "$SESSION_ANALYSED/csf_mask.nii.gz" ]]; do
            sleep 2
        done
    fi

python3 <<PY
import nibabel as nib, numpy as np
from sklearn.decomposition import PCA

ctype = int("$CTYPE")
img = nib.load("$SELFILE")
d = img.get_fdata().astype(np.float32)
X,Y,Z,T = d.shape
flat = d.reshape(-1,T)

def z(x): return (x - x.mean())/(x.std()+1e-6)
regs = []

# -------------------------------------------------------------------------
# 2) Global Signal Regression (GSR)
# -------------------------------------------------------------------------
if ctype == 2:
    regs.append(z(flat.mean(0)))

# -------------------------------------------------------------------------
# 3) tCompCor (top 5% highest variance)
# -------------------------------------------------------------------------
elif ctype == 3:
    std = flat.std(1)
    Nv  = flat.shape[0]
    N   = max(200, int(0.05 * Nv))
    idx = np.argsort(std)[-N:]
    C   = PCA(n_components=5).fit_transform(flat[idx].T)
    for i in range(5): regs.append(z(C[:,i]))

# -------------------------------------------------------------------------
# 4) aCompCor — REQUIRES masks
# -------------------------------------------------------------------------
elif ctype == 4:
    wm  = nib.load("$SESSION_ANALYSED/wm_mask.nii.gz").get_fdata()>0
    csf = nib.load("$SESSION_ANALYSED/csf_mask.nii.gz").get_fdata()>0
    mask = (wm|csf).reshape(-1)
    Xroi = flat[mask].T
    C = PCA(n_components=5).fit_transform(Xroi)
    for i in range(5): regs.append(z(C[:,i]))

# -------------------------------------------------------------------------
# 5) Random CompCor — WORKS WITHOUT MASKS
# -------------------------------------------------------------------------
elif ctype == 5:
    Nv  = flat.shape[0]
    N   = max(200, int(0.05*Nv))
    idx = np.random.choice(Nv, N, replace=False)
    Xroi = flat[idx].T
    C = PCA(n_components=5).fit_transform(Xroi)
    for i in range(5): regs.append(z(C[:,i]))

# -------------------------------------------------------------------------
# 6) Low-variance CompCor — WORKS WITHOUT MASKS
# -------------------------------------------------------------------------
elif ctype == 6:
    std = flat.std(1)
    Nv  = flat.shape[0]
    N   = max(200, int(0.05*Nv))
    idx = np.argsort(std)[:N]
    Xroi = flat[idx].T
    C = PCA(n_components=5).fit_transform(Xroi)
    for i in range(5): regs.append(z(C[:,i]))

# -------------------------------------------------------------------------
# APPLY REGRESSION
# -------------------------------------------------------------------------
if len(regs) > 0:
    R = np.vstack(regs).T
    R = np.hstack([R, np.ones((T,1))])   # intercept
    beta = np.linalg.lstsq(R, flat.T, rcond=None)[0]
    pred = R @ beta
    clean = flat - pred.T
else:
    clean = flat

clean_4d = clean.reshape(X,Y,Z,T)
nib.save(nib.Nifti1Image(clean_4d, img.affine, img.header), "$CONF_O")
print("[OK] Saved confound-regressed → $CONF_O")
PY

    msg_ok "Confound regression complete → $CONF_O"
else
    msg_skip "No confound regression applied."
fi

###############################################################################
# PCA ONLY OPTION
###############################################################################
if [[ "$CTYPE" == "7" ]]; then
    CONF_O="$SELFILE"
fi

###############################################################################
# OPTIONAL PCA DENOISING
###############################################################################
echo -e "${BLUE}Run PCA denoising? (y/n):${RESET}"
read -p "Choice: " RUNPCA

if [[ "$RUNPCA" == "y" ]]; then

    PCA_DIR="${SESSION_ANALYSED}/PCA_${timestamp_now}"
    mkdir -p "$PCA_DIR"
    mkdir -p "${SESSION_ANALYSED}/QC/PCA_grids"

    PCA_OUT="$PCA_DIR/pca_denoised.nii.gz"
    PCA_GRID="${SESSION_ANALYSED}/QC/PCA_grids/pca_grid_${timestamp_now}.png"

python3 <<PY
import numpy as np, nibabel as nib, matplotlib.pyplot as plt
from sklearn.decomposition import PCA

fname = "$CONF_O"
out = "$PCA_OUT"
grid = "$PCA_GRID"

img = nib.load(fname)
d = img.get_fdata().astype(np.float32)
X,Y,Z,T = d.shape
flat = d.reshape(-1,T)

p = PCA(n_components=min(20,T))
C = p.fit_transform(flat.T)

drop=set()
fig,axes=plt.subplots(4,5,figsize=(14,8))
axes=axes.ravel()

def refresh(i):
    ax = axes[i]
    ax.clear()
    t = (C[:,i]-C[:,i].mean())/(C[:,i].std()+1e-6)
    ax.plot(t, 'r' if i in drop else 'k', lw=0.7)
    ax.set_title(f"PC{i+1} {'✗' if i in drop else ''}")
    ax.set_xticks([]); ax.set_yticks([])

for i in range(C.shape[1]):
    refresh(i)

def onclick(ev):
    if ev.inaxes in axes:
        idx = list(axes).index(ev.inaxes)
        if idx in drop: drop.remove(idx)
        else: drop.add(idx)
        refresh(idx)
        fig.canvas.draw_idle()

def onkey(ev):
    if ev.key == "enter":
        plt.savefig(grid, dpi=140)
        print("[OK] Saved PCA grid")
        plt.close()

fig.canvas.mpl_connect("button_press_event", onclick)
fig.canvas.mpl_connect("key_press_event", onkey)

plt.suptitle("Click to DROP components — ENTER to finish")
plt.tight_layout()
plt.show()

for idx in drop:
    C[:,idx] = 0

recon = p.inverse_transform(C).T.reshape(X,Y,Z,T)
nib.save(nib.Nifti1Image(recon.astype(np.float32), img.affine, img.header), out)
print("[OK] PCA denoised saved:", out)
PY

    msg_ok "PCA complete → $PCA_OUT"
    CONF_O="$PCA_OUT"
fi

fi   # end skip step 11

###############################################################################
# ASK USER TO REPEAT STEP 11
###############################################################################
echo
read -p "Repeat Step 11? (y/n): " REP11
if [[ "$REP11" == "y" ]]; then
    exec bash "$0" --continue-from-step 11
fi


###############################################################################
# STEP 12 — SCM (Signal Change Mapping)
###############################################################################
center "STEP 12 — SCM (Signal Change Mapping)"
echo -e "${BOLD}Why?${RESET} Compute baseline, signal, PSC map, and norm_func."

read -p "Skip Step 12? (y/n): " SKIP_12
if [[ "$SKIP_12" != "y" ]]; then

echo -e "${BLUE}Select file for SCM:${RESET}"
mapfile -t FILES < <(find "$SESSION_ANALYSED" -maxdepth 3 -type f -name "*.nii.gz" | sort)
for i in "${!FILES[@]}"; do printf "%3d) %s\n" $((i+1)) "$(basename "${FILES[$i]}")"; done

read -p "Choice: " IDX
SCMFILE="${FILES[$((IDX-1))]}"
SCMNAME=$(basename "$SCMFILE")
msg_ok "Selected → $SCMNAME"

timestamp_scm=$(timestamp)
SCM_DIR="${SESSION_ANALYSED}/SCM_${timestamp_scm}"
mkdir -p "$SCM_DIR"

python3 <<PY
import nibabel as nib, numpy as np, matplotlib.pyplot as plt
img=nib.load("$SCMFILE")
d=img.get_fdata().astype(np.float32)
gs=d.reshape(-1,d.shape[-1]).mean(0)

plt.figure(figsize=(12,4))
plt.plot(gs,'k')
plt.title("Global Mean")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("${SCM_DIR}/global_mean.png",dpi=150)
plt.show()
PY

read -p "Baseline START END volumes: " PB1 PB2
read -p "Signal   START END volumes: " PS1 PS2

python3 <<PY
import nibabel as nib, numpy as np

img=nib.load("$SCMFILE")
d=img.get_fdata()
aff,hdr=img.affine,img.header

b1,b2=int("$PB1"),int("$PB2")
s1,s2=int("$PS1"),int("$PS2")

baseline=d[...,b1:b2+1].mean(-1)
signal=d[...,s1:s2+1].mean(-1)
psc=((signal-baseline)/(baseline+1e-6))*100
norm=(d-baseline[...,None])/(baseline[...,None]+1e-6)

import nibabel as nib
nib.save(nib.Nifti1Image(baseline.astype(np.float32),aff,hdr),f"$SCM_DIR/baseline.nii.gz")
nib.save(nib.Nifti1Image(signal.astype(np.float32),aff,hdr),f"$SCM_DIR/signal.nii.gz")
nib.save(nib.Nifti1Image(psc.astype(np.float32),aff,hdr),f"$SCM_DIR/signal_change_map.nii.gz")
nib.save(nib.Nifti1Image(norm.astype(np.float32),aff,hdr),f"$SCM_DIR/norm_func.nii.gz")
PY

nohup fsleyes "$SCM_DIR/norm_func.nii.gz" "$SCM_DIR/baseline.nii.gz" "$SCM_DIR/signal_change_map.nii.gz" >/dev/null 2>&1 &
msg_ok "SCM saved → $SCM_DIR"

read -p "Repeat Step 12? (y/n): " REP12
if [[ "$REP12" == "y" ]]; then exec "$0" --continue-from-step 12; fi

fi


###############################################################################
# END PART 3
###############################################################################
msg_ok "Part 3/4 complete — ready for QC Summary + CSV LOG."
###############################################################################
# STEP 13 — QC SUMMARY (HTML)
###############################################################################
center "STEP 13 — QC SUMMARY"
echo -e "${BOLD}Why?${RESET} Combine QC metrics, PCA grids, SCM."

QC_HTML="${SESSION_ANALYSED}/QC_summary_$(timestamp).html"

python3 <<PY
import os, glob, datetime

root="$SESSION_ANALYSED"
html=[]

html.append("<html><head><meta charset='UTF-8'>")
html.append("<style>body{background:#111;color:#eee;font-family:Arial;} img{border:2px solid #333;margin:10px;}</style>")
html.append("</head><body>")

html.append("<h1>fUSI QC SUMMARY</h1>")
html.append("<p>Generated: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"</p>")

html.append("<h2>QC Outputs</h2>")
for f in sorted(glob.glob(root+"/QC/*.png")):
    html.append(f"<h3>{os.path.basename(f)}</h3><img src='{f}' width='700'>")

html.append("<h2>SCM Outputs</h2>")
for folder in sorted(glob.glob(root+"/SCM_*")):
    html.append(f"<h3>{os.path.basename(folder)}</h3>")
    for f in sorted(glob.glob(folder+"/*.nii.gz")):
        html.append(f"<p>{os.path.basename(f)}</p>")

html.append("<h2>PCA Grids</h2>")
for grid in sorted(glob.glob(root+"/QC/PCA_grids/*.png")):
    html.append(f"<h3>{os.path.basename(grid)}</h3><img src='{grid}' width='700'>")

html.append("</body></html>")
open("$QC_HTML","w").write("\n".join(html))
print("Saved summary:", "$QC_HTML")
PY

msg_ok "QC Summary saved → $QC_HTML"
###############################################################################
# STEP 14 — PROCESSING LOG (CSV)
###############################################################################
center "STEP 14 — PROCESSING LOG (.CSV)"
echo -e "${BOLD}Why?${RESET} Full reproducibility log."

read -p "Save CSV processing log? (y/n): " SAVELOG
if [[ "$SAVELOG" == "y" ]]; then

    LOG="${SESSION_ANALYSED}/processing_log_$(timestamp).csv"
    echo "Timestamp,File" > "$LOG"

    while IFS= read -r f; do
        rel="${f#$SESSION_ANALYSED/}"
        echo "$(timestamp),$rel" >> "$LOG"
    done < <(find "$SESSION_ANALYSED" -type f | sort)

    msg_ok "Log saved → $LOG"
else
    msg_skip "CSV log skipped."
fi


###############################################################################
# PIPELINE COMPLETE
###############################################################################
center "PIPELINE COMPLETE"
echo -e "${GREEN}All processing steps finished successfully.${RESET}"
echo -e "QC Summary → ${CYAN}$QC_HTML${RESET}"
