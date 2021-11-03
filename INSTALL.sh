wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.4.5_linux64.tar.gz 
tar -zvxf ont-guppy_3.4.5_linux64.tar.gz
mv ont-guppy mop_preprocess/bin/
cd mop_preprocess/bin
ln -s ont-guppy/bin/guppy_* .
ln -s ont-guppy/lib/* .
cd ../../
rm mop_preprocess/bin/ont-guppy/lib/libz.so
cd mop_preprocess/bin/ont-guppy/lib/
ln -s libz.so.1 libz.so
cd ../../../../
rm ont-guppy_3.4.5_linux64.tar.gz
