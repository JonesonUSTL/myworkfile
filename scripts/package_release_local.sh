#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

BUILD_DIR="build-release"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DFEM_ENABLE_EIGEN=OFF -DFEM_ENABLE_PETSC=OFF -DFEM_ENABLE_OPENMP=ON
cmake --build "$BUILD_DIR" -j

mkdir -p dist
OS_NAME="$(uname -s | tr '[:upper:]' '[:lower:]')"
ARCH_NAME="$(uname -m)"
cp "$BUILD_DIR"/fem_solver "dist/fem_solver-${OS_NAME}-${ARCH_NAME}"

tar -czf "fem_solver-${OS_NAME}-${ARCH_NAME}-release.tar.gz" -C dist "fem_solver-${OS_NAME}-${ARCH_NAME}"

echo "Release package created: fem_solver-${OS_NAME}-${ARCH_NAME}-release.tar.gz"
