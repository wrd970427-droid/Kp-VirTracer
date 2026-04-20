#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_ROOT="/root/Kp-VirTracer-test"
LOG_FILE="${TEST_ROOT}/logs/smoke_test.log"

mkdir -p "${TEST_ROOT}/logs"

{
  echo "[INFO] Running environment check"
  bash "${PROJECT_ROOT}/scripts/check_env.sh" || true

  echo "[INFO] Installing package"
  cd "${PROJECT_ROOT}"
  Rscript "${PROJECT_ROOT}/install_package.R"

  echo "[INFO] Running pipeline smoke test"
  Rscript "${PROJECT_ROOT}/run_kpvirtracer.R" \
    --input "${TEST_ROOT}/input" \
    --output "${TEST_ROOT}/output" \
    --threads 4 \
    --config "${TEST_ROOT}/kpvirtracer.test.yaml" \
    --force

  echo "[INFO] Smoke test completed"
} 2>&1 | tee "${LOG_FILE}"
