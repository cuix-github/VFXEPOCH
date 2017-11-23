mkdir $PWD/outputs/anims
mkdir $PWD/outputs/sims

mkdir build && cd build && cmake .. -DVFXEPOCH_EXAMPLES=OFF
make && make install

cd ..

rm -rf build

mkdir build && cd build && cmake .. -DVFXEPOCH_EXAMPLES=ON
make