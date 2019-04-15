# Prepare colvar files
- Compile file2bin.cpp present in src folder with:
    `g++ -O2 ../src/file2bin.cpp -o file2bin.x`
- Convert text files into binary ones and add number of rows and number of columns as first row.
    On each colvar file run:
        `./convert.sh colvar_filename`
- New files, to be used in the input file, will be named as `colvar_filename_new.bin`