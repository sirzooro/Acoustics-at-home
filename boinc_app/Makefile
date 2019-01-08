CAMBALA = ../cambala
ALGLIB = ../alglib-3.14.0/src
BOINC_DIR = ../../boinc
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_ZIP_DIR = $(BOINC_DIR)/zip
FREETYPE_DIR = /usr/include/freetype2

CPPFLAGS = -O3 \
    -std=c++0x -static -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common \
    -DAPP_GRAPHICS -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D NDEBUG \
    -I$(BOINC_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -I$(BOINC_ZIP_DIR) \
    -I$(FREETYPE_DIR) \
    -I$(CAMBALA) \
    -I$(ALGLIB) \
    -L /usr/X11R6/lib \
    -L.

ifeq ($(MinGW32),1)
$(info ===== Compiling MinGW 32-bit app version =====)
CPP = i686-w64-mingw32-g++
else
ifeq ($(MinGW64),1)
$(info ===== Compiling MinGW 64-bit app version =====)
CPP = x86_64-w64-mingw32-g++
else
ifeq ($(M32),1)
$(info ===== Compiling 32-bit app version =====)
CPP = g++ -m32
else
CPP = g++
endif
endif
endif

ifeq ($(SSE2),1)
$(info ===== Compiling SSE2 app version =====)
CPPFLAGS += -msse2
else
ifeq ($(SSE41),1)
$(info ===== Compiling SSE4.1 app version =====)
CPPFLAGS += -msse4.1
else
ifeq ($(AVX),1)
$(info ===== Compiling AVX app version =====)
CPPFLAGS += -mavx -mtune=sandybridge
else
ifeq ($(AVX2),1)
$(info ===== Compiling AVX2 app version =====)
CPPFLAGS += -mavx2 -mfma -mtune=haswell
else
ifeq ($(AVX512),1)
$(info ===== Compiling AVX512 app version =====)
CPPFLAGS += -march=skylake-avx512
ifeq ($(MinGW64),1)
# MinGW needs workaround for "invalid register for .seh_savexmm" bug
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=65782
CPPFLAGS += \
	-ffixed-xmm16 -ffixed-xmm17 -ffixed-xmm18 -ffixed-xmm19 -ffixed-xmm20 -ffixed-xmm21 -ffixed-xmm22 -ffixed-xmm23 \
	-ffixed-xmm24 -ffixed-xmm25 -ffixed-xmm26 -ffixed-xmm27 -ffixed-xmm28 -ffixed-xmm29 -ffixed-xmm30 -ffixed-xmm31
endif
endif
endif
endif
endif
endif

ifeq ($(NOWARN),1)
CPPFLAGS += -w
endif

release: cambala_boinc_app

libstdc++.a:
	ln -s `${CPP} -print-file-name=libstdc++.a`

clean: distclean

distclean:
	/bin/rm -f $(PROGS) *.o libstdc++.a cambala_boinc_app cambala_boinc_app.exe

cambala_boinc_app: cambala_boinc_app.o all.o \
	libstdc++.a $(BOINC_API_DIR)/libboinc_api.a $(BOINC_LIB_DIR)/libboinc.a
	$(CPP) $(CPPFLAGS) cambala_boinc_app.o all.o -o cambala_boinc_app \
	libstdc++.a -pthread $(BOINC_API_DIR)/libboinc_api.a \
	$(BOINC_LIB_DIR)/libboinc.a
	
all.o: all.cpp $(ALGLIB)/specialfunctions.cpp $(ALGLIB)/linalg.cpp $(ALGLIB)/ap.cpp $(ALGLIB)/alglibmisc.cpp $(ALGLIB)/alglibinternal.cpp $(CAMBALA)/sequential.cpp
	${CPP} ${CPPFLAGS} all.cpp -c

cambala_boinc_app.o: cambala_boinc_app.cpp
	${CPP} ${CPPFLAGS} cambala_boinc_app.cpp -c
