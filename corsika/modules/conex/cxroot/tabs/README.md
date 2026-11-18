# corsika-data

New place to keep all the CORSIKA7, CORSIKA8, CONEX, CRMC etc. large
data files.

The project contains cmake build files. It can be included in our project to access corsika-data by adding
`add_subdirectory (corsika-data)` to your CMakeLists.txt. 
This creates the cmake taget "CorsikaData" and it builds the libCorsikaData.a archive. 


Update, July 2020, RU
---------------------
The larger data files are now bzip2-ed There is a helper lib 'readLib'
whith cmake file that can be used for convenience to unpack and read
those files directly. Fortran wrappers are provided, too. 
