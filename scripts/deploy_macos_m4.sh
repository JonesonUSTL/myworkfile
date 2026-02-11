#!/usr/bin/env bash
set -euo pipefail

# macOS Apple Silicon(M4/ARM64) 一键部署脚本。
# 功能：检查 Homebrew -> 安装依赖 -> CMake 配置与编译 -> 运行冒烟用例。

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

echo "[1/5] 检查平台..."
uname -s | grep -q "Darwin" || { echo "仅支持 macOS"; exit 1; }
uname -m | grep -Eq "arm64|aarch64" || echo "警告：当前非 ARM64，脚本仍继续。"

echo "[2/5] 检查 Homebrew..."
if ! command -v brew >/dev/null 2>&1; then
  echo "未检测到 Homebrew，请先安装：https://brew.sh"
  exit 1
fi

echo "[3/5] 安装依赖..."
brew install cmake eigen || true

echo "[4/5] 配置与编译..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFEM_ENABLE_EIGEN=ON -DFEM_ENABLE_PETSC=OFF -DFEM_ENABLE_OPENMP=ON
cmake --build build -j

echo "[5/5] 运行冒烟测试..."
./build/fem_solver examples/linear_truss.inp output/linear_macos.vtk 1.0 dense

echo "部署完成。可执行文件：build/fem_solver"
