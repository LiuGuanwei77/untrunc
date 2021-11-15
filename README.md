Untrunc
=======

This is a fork of Federico Ponchio's repository: https://github.com/ponchio/untrunc. Untrunc is a tool to restore a damaged (truncated) mp4, m4v, mov, 3gp video and m4a audio.

Changes to the original repository:

* Change the dependency from Libav to Ffmpeg(https://www.ffmpeg.org)
* Add support for Windows


## Installing on Ubuntu

1. Install CMake

2. Clone Ffmpeg's repository
   ```bash
   git clone https://git.ffmpeg.org/ffmpeg.git ffmpeg
   ```

3. Build Ffmpeg
   ```bash
   cd ffmpeg
   ./configure --enable-shared
   make
   ```
4. Generate Untrunc's Makefile
   ```bash
   cmake -G "Unix Makefiles" -B "unix_build"
   ```

5. Edit unix_build/CMakeCache.txt, set the value of variable AVCODEC_LIBRARY, AVFORMAT_LIBRARY, AVUTIL_LIBRARY, FFMPEG_INCLUDE_DIR, FFMPEG_SOURCE_DIR

6. Build Untrunc
   ```bash
   cd unix_build
   make
   ```

## Installing on Windows

1. Install CMake

   https://cmake.org/download/

2. Clone Ffmpeg's repository
   ```bash
   git clone https://git.ffmpeg.org/ffmpeg.git ffmpeg
   ```

3. Build Ffmpeg

   https://trac.ffmpeg.org/wiki/CompilationGuide/MSVC

4. Generate Untrunc's Makefile
   Open CMake GUI, set the source directory and the build directory.
   Push "Configure" button.

5. Set the value of variable AVCODEC_LIBRARY, AVFORMAT_LIBRARY, AVUTIL_LIBRARY, FFMPEG_INCLUDE_DIR, FFMPEG_SOURCE_DIR. Push "Generate" button.

5. Open untrunc.sln with Visual Studio and build


## Using

You need both the broken video and an example working video (ideally from the same camera, if not the chances to fix it are slim).

Run this command in the folder where you have unzipped and compiled Untrunc but replace the `/path/to/...` bits with your 2 video files:
   ```bash
   ./untrunc /path/to/working-video.m4v /path/to/broken-video.m4v
   ```
Then it should churn away and hopefully produce a playable file called `broken-video_fixed.m4v`.

That's it you're done!

