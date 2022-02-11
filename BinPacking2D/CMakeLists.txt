cmake_minimum_required(VERSION 3.20)
project(BinPacking2D)

set(CMAKE_CXX_STANDARD 11)

include_directories(3rdParty/brkgaAPI-src/src 3rdParty/RectangleBinPack/include)
add_executable(BinPacking2D main.cpp
        3rdParty/RectangleBinPack/GuillotineBinPack.cpp
        3rdParty/RectangleBinPack/MaxRectsBinPack.cpp
        3rdParty/RectangleBinPack/Rect.cpp
        3rdParty/RectangleBinPack/ShelfBinPack.cpp
        3rdParty/RectangleBinPack/ShelfNextFitBinPack.cpp
        3rdParty/RectangleBinPack/SkylineBinPack.cpp
)