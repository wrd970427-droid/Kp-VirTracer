#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${PROJECT_ROOT}/environment.yml"
ENV_NAME="kpvirtracer"

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda not found. Please install Miniconda/Anaconda first."
  exit 1
fi

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "[ERROR] Missing ${ENV_FILE}"
  exit 1
fi

env_exists() {
  conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}" \
    || [[ -d "$(conda info --base)/envs/${ENV_NAME}" ]]
}

echo "[INFO] Checking conda environment: ${ENV_NAME}"
if env_exists; then
  echo "[INFO] Environment '${ENV_NAME}' already exists; skipping create/update."
  echo "[INFO] To recreate from scratch: conda env remove -n ${ENV_NAME} -y && bash \"${PROJECT_ROOT}/install_env.sh\""
  echo "[INFO] To sync packages with environment.yml yourself: conda env update -n ${ENV_NAME} -f \"${ENV_FILE}\" --prune"
else
  echo "[INFO] Creating environment '${ENV_NAME}' from environment.yml ..."
  conda env create -f "${ENV_FILE}"
fi

cat <<'EOF'
[INFO] Environment setup complete.
Next steps:
  1) conda activate kpvirtracer
  2) Rscript install_package.R
  3) bash scripts/check_env.sh
  4) Rscript run_kpvirtracer.R --help
EOF
