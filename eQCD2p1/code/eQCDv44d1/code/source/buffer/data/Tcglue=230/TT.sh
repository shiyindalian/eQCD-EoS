#!/bin/bash
for a in {0..80} 
do
cp -r Tem$a/buffer data/mu$a
done
#rm -rf Tem*
