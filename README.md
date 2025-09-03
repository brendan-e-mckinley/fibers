# fibers


compilation: `g++ -arch x86_64h $(python3-config --cflags) $(python3-config --ldflags) new\ code/c_fibers_obj_clean.cpp \
-I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/ \
-shared -undefined dynamic_lookup -std=c++17 \
-o c_fibers_obj_clean.so`

c++ -O3 -shared new\ code/c_fibers_obj_clean.cpp -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/ -o c_fibers_obj_clean.so -std=c++14 -fPIC -DDOUBLE_PRECISION `python3 -m pybind11 --includes` 