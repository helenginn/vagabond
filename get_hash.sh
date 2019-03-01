#!/bin/bash

if [ $# -eq 0 ]; then
	find libsrc libgui -name '*.cpp' -o -name '*.[c:h]' > get_hash.dep
	exit
fi

vagabond_version=`grep 'VAGABOND_VERSION_NUMBER' $1 | cut -d' ' -f3`
vagabond_commit="notcompiledfromgit"

command -v git > /dev/null 2>&1
if [ $? -eq 0 ]; then
	vagabond_commit="."`git rev-parse HEAD`
fi

sed "s/_version_commit_/${vagabond_version}${vagabond_commit}/g" $1 > $2


