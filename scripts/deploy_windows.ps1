# Windows PowerShell 一键部署脚本
# 功能：CMake 配置与构建，并执行线性桁架冒烟算例。

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

Write-Host "[1/3] 配置项目..."
cmake -S . -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DFEM_ENABLE_EIGEN=ON -DFEM_ENABLE_PETSC=OFF -DFEM_ENABLE_OPENMP=ON

Write-Host "[2/3] 编译..."
cmake --build build -j

Write-Host "[3/3] 运行冒烟测试..."
.\build\fem_solver.exe examples\linear_truss.inp output\linear_windows.vtk 1.0 dense

Write-Host "部署完成。可执行文件：build\\fem_solver.exe"
