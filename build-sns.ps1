param([string]$BuildDirectory = (Join-Path $PSScriptRoot 'build\sns'))
$ErrorActionPreference = 'Stop'
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
$vs = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vs) { throw 'Visual Studio C++ toolchain not found' }
$msvc = Get-ChildItem -LiteralPath (Join-Path $vs 'VC\Tools\MSVC') -Directory | Sort-Object Name -Descending | Select-Object -First 1
$kits = Join-Path ${env:ProgramFiles(x86)} 'Windows Kits\10'
$sdk = Get-ChildItem -LiteralPath (Join-Path $kits 'Include') -Directory | Sort-Object Name -Descending | Select-Object -First 1
$oldInclude = $env:INCLUDE
$oldLib = $env:LIB
try {
    $env:INCLUDE = @((Join-Path $msvc.FullName 'include'),(Join-Path $sdk.FullName 'ucrt'),(Join-Path $sdk.FullName 'shared'),(Join-Path $sdk.FullName 'um')) -join ';'
    $env:LIB = @((Join-Path $msvc.FullName 'lib\x64'),(Join-Path $kits "Lib\$($sdk.Name)\ucrt\x64"),(Join-Path $kits "Lib\$($sdk.Name)\um\x64")) -join ';'
    New-Item -ItemType Directory -Path $BuildDirectory -Force | Out-Null
    $destination = (Resolve-Path -LiteralPath $BuildDirectory).Path
    $compiler = Join-Path $msvc.FullName 'bin\Hostx64\x64\cl.exe'
    $source = Join-Path $PSScriptRoot 'Keldysh\sns'
    $common = @('/nologo','/std:c++17','/EHsc','/O2','/W4','/utf-8',"/Fo$destination\",(Join-Path $source 'spectral.cpp'),(Join-Path $source 'kinetic.cpp'),(Join-Path $source 'current.cpp'),(Join-Path $source 'runner.cpp'))
    & $compiler @common (Join-Path $source 'tests.cpp') "/Fe$destination\sns_tests.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Test build failed' }
    & $compiler @common (Join-Path $source 'cli.cpp') "/Fe$destination\sns_cli.exe"
    if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
    & $compiler /nologo /std:c++17 /EHsc /O2 /W4 /utf-8 "/Fo$destination\" (Join-Path $source 'anderson_safeguard_tests.cpp') "/Fe$destination\sns_anderson_safeguard_tests.exe"
    if ($LASTEXITCODE -ne 0) { throw 'Anderson safeguard build failed' }
    & (Join-Path $destination 'sns_anderson_safeguard_tests.exe')
    if ($LASTEXITCODE -ne 0) { throw 'Anderson safeguard tests failed' }
    & (Join-Path $destination 'sns_tests.exe')
    if ($LASTEXITCODE -ne 0) { throw 'Tests failed' }
} finally { $env:INCLUDE = $oldInclude; $env:LIB = $oldLib }
