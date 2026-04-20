#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="kpvirtracer"
TOOLS=(kleborate mob_recon blastn makeblastdb fastANI prodigal Rscript)

echo "== Kp-VirTracer Environment Check =="

if command -v conda >/dev/null 2>&1; then
  echo "[OK] conda found: $(command -v conda)"
else
  echo "[FAIL] conda not found"
  exit 1
fi

active_env="${CONDA_DEFAULT_ENV:-<none>}"
if [[ "${active_env}" == "${ENV_NAME}" ]]; then
  echo "[OK] active conda env: ${active_env}"
else
  echo "[WARN] active conda env is '${active_env}', expected '${ENV_NAME}'"
fi
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}" || [[ -d "${HOME}/miniconda3/envs/${ENV_NAME}" ]]; then
  echo "[OK] conda env '${ENV_NAME}' exists"
else
  echo "[WARN] conda env '${ENV_NAME}' not found"
fi

printf "\n%-15s %-8s %s\n" "Tool" "Status" "Path"
printf "%-15s %-8s %s\n" "----" "------" "----"
for t in "${TOOLS[@]}"; do
  if command -v "${t}" >/dev/null 2>&1; then
    printf "%-15s %-8s %s\n" "${t}" "OK" "$(command -v "${t}")"
  else
    printf "%-15s %-8s %s\n" "${t}" "MISSING" "-"
  fi
done
