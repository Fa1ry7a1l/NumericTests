$ErrorActionPreference = "Stop"

$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$mingwBin = "C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin"
$gitUsrBin = "C:\PROGRA~1\Git\usr\bin"
$gitBin = "C:\Program Files\Git\bin"
$toolsDir = Join-Path $root "tools"
$msysRoot = "/" + ($root.Substring(0, 1).ToLower()) + $root.Substring(2).Replace("\", "/")
$msysTools = "$msysRoot/tools"
$msysNtlSrc = "$msysRoot/third_party/ntl/src"
$msysM4Dir = "$msysRoot/third_party/msys-tools/usr/bin"

New-Item -ItemType Directory -Force -Path $toolsDir | Out-Null
Copy-Item -Force -LiteralPath (Join-Path $mingwBin "mingw32-make.exe") -Destination (Join-Path $toolsDir "make.exe")

$env:PATH = "$toolsDir;$mingwBin;$gitUsrBin;$gitBin;$env:PATH"
$sh = Join-Path $gitUsrBin "sh.exe"

$bashPath = "'${msysTools}:${msysM4Dir}:/c/PROGRA~1/JETBRA~1/CLION2~1.3/bin/mingw/bin:/c/PROGRA~1/Git/usr/bin:/usr/bin:/bin'"

if (-not (Test-Path (Join-Path $root "third_party\gmp-install\lib\libgmp.a"))) {
    & (Join-Path $root "scripts\build_gmp.ps1")
    New-Item -ItemType Directory -Force -Path $toolsDir | Out-Null
    Copy-Item -Force -LiteralPath (Join-Path $mingwBin "mingw32-make.exe") -Destination (Join-Path $toolsDir "make.exe")
}

& $sh -lc "export SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export CONFIG_SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export PATH=${bashPath}:`$PATH; cd '$msysNtlSrc'; ./configure PREFIX=../../ntl-install GMP_PREFIX=../../gmp-install NTL_GMP_LIP=on SHARED=off NTL_THREADS=off NTL_EXCEPTIONS=on CXX=g++"
& $sh -lc "export SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export CONFIG_SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export PATH=${bashPath}:`$PATH; cd '$msysNtlSrc'; make -j4"
& $sh -lc "export SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export CONFIG_SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export PATH=${bashPath}:`$PATH; cd '$msysNtlSrc'; make install"

Remove-Item -Force -LiteralPath (Join-Path $toolsDir "make.exe") -ErrorAction SilentlyContinue
if ((Get-ChildItem -Force $toolsDir -ErrorAction SilentlyContinue | Measure-Object).Count -eq 0) {
    Remove-Item -Force -LiteralPath $toolsDir
}
