# Vendored: mindis/hll

Source: https://github.com/mindis/hll (commit e061f27, "Compile both C and C++ files with C++ compiler.")
Vendored on: 2026-05-06
Files copied: hll.cpp, hll.h, kthread.c, kthread.h, logutil.h, sseutil.h, LICENSE.

## Why vendored, not a submodule
- hll has no recent upstream activity, so a submodule buys little.
- Vendoring keeps the hammock_claude build offline-friendly.

## Updating
If you need to refresh: clone https://github.com/mindis/hll, copy the same
file set in, and update the commit reference above. CMakeLists.txt in this
directory is local to hammock_claude (not from upstream) and should be
preserved across updates.
