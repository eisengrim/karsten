#!/bin/bash

GIT_VERSION=$(git describe --always --dirty="X" --tags)
GIT_BRANCH=$(git symbolic-ref -q --short HEAD)
if [ -z "${GIT_BRANCH}" ]
then
	GIT_VERSION_STRING=${GIT_VERSION}
else
	GIT_VERSION_STRING="${GIT_VERSION}-${GIT_BRANCH}"
fi

GIT_VERSION_MESH=$(git log --pretty=format:%H -1 mesh.h mesh.c)
GIT_VERSION_PARTICLES=$(git log --pretty=format:%H -1 particles.h particles.c)
 
GIT_VERSION_H="// Version Information from git \n\
#ifndef F_VERSION_H\n\
#define F_VERSION_H\n\
#define GIT_VERSION_STRING \"${GIT_VERSION_STRING}\"\n\
#define GIT_VERSION_MESH \"${GIT_VERSION_MESH}\"\n\
#define GIT_VERSION_PARTICLES \"${GIT_VERSION_PARTICLES}\"\n\
#endif"

if [ ! -f version.h ]
then 
	touch version.h
fi

echo -e "${GIT_VERSION_H}" | cmp -s version.h -
R=$?
if [ "$R" -gt "0" ]
then
	echo -e "${GIT_VERSION_H}" > version.h
	gcc -x c version.h -c -o version.o
fi	

