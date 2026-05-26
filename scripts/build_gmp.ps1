$ErrorActionPreference = "Stop"

$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$mingwBin = "C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin"
$gitUsrBin = "C:\PROGRA~1\Git\usr\bin"
$toolsDir = Join-Path $root "tools"
$thirdParty = Join-Path $root "third_party"
$msysTools = Join-Path $thirdParty "msys-tools"
$gmpArchive = Join-Path $thirdParty "gmp-6.3.0.tar.xz"
$m4Archive = Join-Path $msysTools "m4.pkg.tar.zst"

New-Item -ItemType Directory -Force -Path $toolsDir, $thirdParty, $msysTools | Out-Null
Copy-Item -Force -LiteralPath (Join-Path $mingwBin "mingw32-make.exe") -Destination (Join-Path $toolsDir "make.exe")

if (-not (Test-Path $gmpArchive)) {
    curl.exe -L --retry 3 --fail https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz -o $gmpArchive
}

if (-not (Test-Path (Join-Path $msysTools "usr\bin\m4.exe"))) {
    curl.exe -L --retry 3 --fail https://repo.msys2.org/msys/x86_64/m4-1.4.21-1-x86_64.pkg.tar.zst -o $m4Archive
    tar -xf $m4Archive -C $msysTools
}

if (-not (Test-Path (Join-Path $thirdParty "gmp-6.3.0"))) {
    tar -xf $gmpArchive -C $thirdParty
}

$msysRoot = "/" + ($root.Substring(0, 1).ToLower()) + $root.Substring(2).Replace("\", "/")
$msysToolsDir = "$msysRoot/tools"
$msysM4Dir = "$msysRoot/third_party/msys-tools/usr/bin"
$msysGmpSrc = "$msysRoot/third_party/gmp-6.3.0"
$msysGmpInstall = "$msysRoot/third_party/gmp-install"
$bashPath = "'${msysToolsDir}:${msysM4Dir}:/c/PROGRA~1/JETBRA~1/CLION2~1.3/bin/mingw/bin:/c/PROGRA~1/Git/usr/bin:/usr/bin:/bin'"

$env:PATH = "$toolsDir;$(Join-Path $msysTools "usr\bin");$mingwBin;$gitUsrBin;$env:PATH"
$sh = Join-Path $gitUsrBin "sh.exe"

& $sh -lc "export SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export CONFIG_SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export PATH=${bashPath}:`$PATH; cd '$msysGmpSrc'; make distclean >/dev/null 2>&1 || true; rm -rf '$msysGmpInstall'"
& $sh -lc "export SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export CONFIG_SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export PATH=${bashPath}:`$PATH; cd '$msysGmpSrc'; ./configure --prefix='$msysGmpInstall' --enable-static --disable-shared --with-pic --enable-cxx ABI=64 CC=gcc CXX=g++ NM=nm LD=ld AR=ar RANLIB=ranlib"
& $sh -lc "export SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export CONFIG_SHELL=/c/PROGRA~1/Git/usr/bin/sh.exe; export PATH=${bashPath}:`$PATH; cd '$msysGmpSrc'; make -j4"

New-Item -ItemType Directory -Force -Path (Join-Path $thirdParty "gmp-install\include"), (Join-Path $thirdParty "gmp-install\lib") | Out-Null
Copy-Item -Force (Join-Path $thirdParty "gmp-6.3.0\gmp.h") (Join-Path $thirdParty "gmp-install\include\gmp.h")
Copy-Item -Force (Join-Path $thirdParty "gmp-6.3.0\gmpxx.h") (Join-Path $thirdParty "gmp-install\include\gmpxx.h")
Copy-Item -Force (Join-Path $thirdParty "gmp-6.3.0\.libs\libgmp.a") (Join-Path $thirdParty "gmp-install\lib\libgmp.a")
Copy-Item -Force (Join-Path $thirdParty "gmp-6.3.0\.libs\libgmpxx.a") (Join-Path $thirdParty "gmp-install\lib\libgmpxx.a")

Remove-Item -Force -LiteralPath (Join-Path $toolsDir "make.exe") -ErrorAction SilentlyContinue
if ((Get-ChildItem -Force $toolsDir -ErrorAction SilentlyContinue | Measure-Object).Count -eq 0) {
    Remove-Item -Force -LiteralPath $toolsDir
}
