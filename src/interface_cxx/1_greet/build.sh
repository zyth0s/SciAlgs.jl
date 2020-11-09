
rm -Rf build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/Users/daniel/.julia/artifacts/6017255205dc4fbf4d962903a855a0c631f092dc ..
cmake --build . --config Release
