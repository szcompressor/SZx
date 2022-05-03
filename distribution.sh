#!/bin/bash

packageName=szx-1.0
minor=0
revision=0
versiontype=stable

tarName=${packageName}.${minor}.${revision}
newPackName=${packageName}.${minor}.${revision}
rm -rf ${newPackName}
make dist

tar -xzvf ${packageName}.tar.gz
cp -rf example/README.md example/Makefile.bk ${packageName}/example/
cp CMakeLists.txt ${packageName}
cp example/CMakeLists.txt ${packageName}/example/
cp szx/CMakeLists.txt ${packageName}/szx
cp szx.pc.in ${packageName}
cp copyright-and-BSD-license.txt ${packageName}

mv ${packageName} ${newPackName}

rm ${packageName}.tar.gz
tar -czvf ${tarName}.tar.gz ${newPackName}


