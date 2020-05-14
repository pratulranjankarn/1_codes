#!/bin/bash

for i in ~/COMB_EV/Event_wav_final/2*
do
cd "${i}"
echo ${i}
rm -f out*.txt
cd ..
done
