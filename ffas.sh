#!/bin/bash
#environment
export FFAS=/Folder/of/FFAS/installation
export PATH=$PATH:/Folder/of/FFAS/installation/soft
export MAX_PPN=12
export SCRATCH_DISK= /Folder/for/temporary/files

cd /path/to/folder/of/this/script


#preparing profiles
ff $2 $2pr.2 $2mu.2 $2ff.2
binff $2ff.2 $2fb.2

#comparing them
ffasn 2OMZ.ff $2fb.2 -a > $3

mv a.txt $4
rm -f $2pr.2
rm -f $2mu.2
rm -f $2ff.2
rm -f $2fb.2

cd -