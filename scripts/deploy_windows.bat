@echo off
setlocal

REM Windows CMD 一键部署脚本。
cd /d %~dp0\..

echo [1/3] Configure...
cmake -S . -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DFEM_ENABLE_EIGEN=ON -DFEM_ENABLE_PETSC=OFF -DFEM_ENABLE_OPENMP=ON || exit /b 1

echo [2/3] Build...
cmake --build build -j || exit /b 1

echo [3/3] Smoke test...
build\fem_solver.exe examples\linear_truss.inp output\linear_windows.vtk 1.0 dense || exit /b 1

echo Done. Binary: build\fem_solver.exe
endlocal
