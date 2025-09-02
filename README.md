# fibers


compilation: `gcc $(python3-config --cflags) $(python3-config --ldflags) new\ code/c_fibers_obj_clean.cpp -o c_fibers_obj_clean.so -std=c++17 -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/ -shared -undefined dynamic_lookup`