LIBS = `root-config --libs` `root-config --auxlibs` -lTMVA
CFLAGS = `root-config --cflags`
CC = g++

framework:  MyDict.cxx
	$(CC) framework.cc MyDict.cxx -o run_framework $(LIBS) $(CFLAGS)

MyDict.cxx:  Linkdef.h
	rootcint -f $@ -c -p $^

%_cint.cxx: %.h %_Linkdef.h
	rootcint -f $@ -c $*.h $*_Linkdef.h
