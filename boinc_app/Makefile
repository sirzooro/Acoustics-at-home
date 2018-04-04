CAMBALA = ..
ALGLIB = ../../alglib
BOINC_DIR = ../../boinc
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_ZIP_DIR = $(BOINC_DIR)/zip
FREETYPE_DIR = /usr/include/freetype2

CPP = g++
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

release: cambala_boinc

libstdc++.a:
	ln -s `g++ -print-file-name=libstdc++.a`

ttfont.cpp:
	ln -s ../../api/ttfont.cpp .

clean: distclean

distclean:
	/bin/rm -f $(PROGS) *.o libstdc++.a cambala_boinc

cambala_boinc: cambala_boinc.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o cambala_sequential.o \
	libstdc++.a $(BOINC_API_DIR)/libboinc_api.a $(BOINC_LIB_DIR)/libboinc.a
	$(CPP) $(CPPFLAGS) cambala_boinc.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o cambala_sequential.o -o cambala_boinc \
	libstdc++.a -pthread $(BOINC_API_DIR)/libboinc_api.a \
	$(BOINC_LIB_DIR)/libboinc.a
	
specialfunctions.o: $(ALGLIB)/specialfunctions.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/specialfunctions.cpp -c

linalg.o: $(ALGLIB)/linalg.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/linalg.cpp -c

ap.o: $(ALGLIB)/ap.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/ap.cpp -c

alglibmisc.o: $(ALGLIB)/alglibmisc.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/alglibmisc.cpp -c

alglibinternal.o: $(ALGLIB)/alglibinternal.cpp 
	${CPP} ${CPPFLAGS} $(ALGLIB)/alglibinternal.cpp -c

cambala_sequential.o: $(CAMBALA)/cambala_sequential.cpp
	${CPP} ${CPPFLAGS} $(CAMBALA)/cambala_sequential.cpp -c

cambala_boinc.o: cambala_boinc.cpp
	${CPP} ${CPPFLAGS} cambala_boinc.cpp -c