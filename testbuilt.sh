mkdir -p build/temp.macosx-10.9-x86_64-2.7/src/nest
mkdir -p build/temp.macosx-10.9-x86_64-2.7/src/util

clang++ -stdlib=libstdc++ -mmacosx-version-min=10.6 -fno-strict-aliasing -fno-common -dynamic -I/usr/local/include -I/usr/local/opt/sqlite/include -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -I/usr/local/Cellar/eigen/3.2.1/include/eigen3 -I/usr/include/c++/4.2.1 -Isrc -I/usr/local/Cellar/python/2.7.6/Frameworks/Python.framework/Versions/2.7/include/python2.7 -c src/nest/NEST_wrap.cpp -o build/temp.macosx-10.9-x86_64-2.7/src/nest/NEST_wrap.o

clang++ -stdlib=libstdc++ -mmacosx-version-min=10.6 -bundle -undefined dynamic_lookup -L/usr/local/lib -L/usr/local/opt/sqlite/lib build/temp.macosx-10.9-x86_64-2.7/src/nest/NEST_wrap.o -o /Users/sheffler/Dropbox/project/scheme/scheme/NEST_wrap.so

clang++ -stdlib=libstdc++ -mmacosx-version-min=10.6 -fno-strict-aliasing -fno-common -dynamic -I/usr/local/include -I/usr/local/opt/sqlite/include -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -I/usr/local/Cellar/eigen/3.2.1/include/eigen3 -I/usr/include/c++/4.2.1 -Isrc -I/usr/local/Cellar/python/2.7.6/Frameworks/Python.framework/Versions/2.7/include/python2.7 -c src/util/util_tests.cpp -o build/temp.macosx-10.9-x86_64-2.7/src/util/util_tests.o

clang++ -stdlib=libstdc++ -mmacosx-version-min=10.6 -bundle -undefined dynamic_lookup -L/usr/local/lib -L/usr/local/opt/sqlite/lib build/temp.macosx-10.9-x86_64-2.7/src/util/util_tests.o -o /Users/sheffler/Dropbox/project/scheme/scheme/util_tests.so
