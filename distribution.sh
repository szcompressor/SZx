#!/bin/bash

packageName=sz-2.1
minor=12
revision=0
versiontype=stable

#tarName=${packageName}.${minor}-${versiontype}
tarName=${packageName}.${minor}.${revision}
#newPackName=${packageName}.${minor}.${revision}-${versiontype}
newPackName=${packageName}.${minor}.${revision}
#rm -rf ${packageName}
rm -rf ${newPackName}
make dist

tar -xzvf ${packageName}.tar.gz
cp COPYRIGHT.txt ${packageName}/
cp -rf example/README.md example/Makefile.bk ${packageName}/example/
cp -rf doc ${packageName}/
rm -rf ${packageName}/sz/*.mod
rm -rf example/testdata/x86/*.sz example/testdata/x86/*.out example/testdata/x86/*.h
cp -r example/testdata ${packageName}/example/
cp example/sz.config ${packageName}/example/
cp example/sz_int.config ${packageName}/example/
cp example/test.sh ${packageName}/example/
cp example/test_int.sh ${packageName}/example/
cp example/testfloat_CompDecomp*.c ${packageName}/example/
cp example/testdouble_CompDecomp*.c ${packageName}/example/
cp example/compile-CompDecomp.sh ${packageName}/example/
cp -r test ${packageName}/
cd hdf5-filter/H5Z-SZ
make clean
cd test
make clean
rm *.h5
cd ../../..
cp -r hdf5-filter ${packageName}/
cp -r adiosReader ${packageName}/
#copy cmake files
cp CMakeLists.txt ${packageName}/
cp sz/CMakeLists.txt ${packageName}/sz
cp zlib/CMakeLists.txt ${packageName}/zlib
cp zstd/CMakeLists.txt ${packageName}/zstd
cp example/CMakeLists.txt ${packageName}/example
cp -r cmake ${packageName}/

cd NetCDFReader
make clean -f Makefile.bk
cd ..
cp -r NetCDFReader ${packageName}/

mv ${packageName} ${newPackName}

rm ${packageName}.tar.gz
tar -czvf ${tarName}.tar.gz ${newPackName}


