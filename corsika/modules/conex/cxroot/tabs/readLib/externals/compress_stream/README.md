# compress_stream

Header only c++ library for read and write compress file with std like stream

## Credits

Code comes from [nice-core project](https://github.com/cvjena/nice-core), I just merge declaration and implementation of class in one header file, and add test.

## How to use

### Instalation

#### Git submodule

Add compress_stream as submodule :

```
git submodule add https://github.com/natir/compress_stream.git
```

#### Hard copy

Download include file :

```
wget https://raw.githubusercontent.com/natir/compress_stream/master/include/gzstream.hpp
wget https://raw.githubusercontent.com/natir/compress_stream/master/include/bzstream.hpp
```

### Compile with cmake

```cmake
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

find_package(BZip2 REQUIRED)
include_directories(${BZIP2_INCLUDE_DIR})


# If you use git submodule uncomment this line
# include_directories(compress_stream/include/)
# else
# include_directories(path/to/header/file/)

file(
    GLOB_RECURSE
    src_files
    src/*
    )

add_executable(${executable_name} ${src_files})
target_link_libraries(${executable_name} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})
```

### Compile with manual build system

```
g++ $(pkg-config --cflags --libs zlib) -lbz2 dev/compress_stream/test/src/main.cpp -I path/to/header/file/
```


## Documentation

Use (i|o)gzstream and (i|o)bzstream exactly like classic (i|o)stream.

