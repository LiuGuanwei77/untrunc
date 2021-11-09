
cmake_minimum_required(VERSION 3.12)

# Only generate Debug and Release configuration types.
set(CMAKE_CONFIGURATION_TYPES Debug Release)

# Project name.
# TODO: Change this line to match your project name when you copy this file.
project(untrunc VERSION 1.0.0 LANGUAGES CXX)

set(UNTRUNC_HDRS
    src/AP_AtomDefinitions.h
    src/atom.h
    src/codec.h
    src/codecstats.h
    src/file.h
    src/log.h
    src/mp4.h
    src/track.h
)
set(UNTRUNC_SRCS
    src/atom.cpp
    src/codec.cpp
    src/codec_alac.cpp
    src/codec_apch.cpp
    src/codec_avc1.cpp
    src/codec_camm.cpp
    src/codec_fdsc.cpp
    src/codec_gpmd.cpp
    src/codec_hev1.cpp
    src/codec_mbex.cpp
    src/codec_mijd.cpp
    src/codec_mp4a.cpp
    src/codec_mp4v.cpp
    src/codec_pcm.cpp
    src/codec_rtp.cpp
    src/codec_text.cpp
    src/codec_tmcd.cpp
    src/codec_unknown.cpp
    src/codecstats.cpp
    src/file.cpp
    src/log.cpp
    src/main.cpp
    src/mp4.cpp
    src/track.cpp
) 

if(MSVC)
    add_compile_definitions(NOMINMAX) 
endif()

add_executable(untrunc ${UNTRUNC_SRCS} ${UNTRUNC_HDRS})

find_path(FFMPEG_INCLUDE_DIR libavcodec/avcodec.h)
find_library(AVCODEC_LIBRARY avcodec)

find_library(AVFORMAT_LIBRARY avformat)

find_library(AVUTIL_LIBRARY avutil)

target_include_directories(untrunc PUBLIC ${FFMPEG_INCLUDE_DIR})
#target_link_directories(untrunc PUBLIC ${FFMPEG_LIBDIR})
target_link_libraries(untrunc ${AVCODEC_LIBRARY} ${AVFORMAT_LIBRARY} ${AVUTIL_LIBRARY})