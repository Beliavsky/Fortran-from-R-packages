# Backend directory

`backend/bin` is used for the Windows DLL and `backend/lib` for Unix shared libraries.

Generated CMake build trees and installed static libraries are not included in release archives. Run `scripts/build_backend.bat` or `scripts/build_backend.sh` to recreate them from the bundled OSQP sources.
