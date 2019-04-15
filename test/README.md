# Prepare colvar files
- Compile file2bin.cpp present in src folder with:
    `g++ -O2 ../src/file2bin.cpp -o file2bin.x`
- Convert text files into binary ones and add number of rows and number of columns as first row.
    On each colvar file run:
        `./convert.sh colvar_filename`
- New files, to be used in the input file, will be named as `colvar_filename_new.bin`

***
**NOTE**: For storage reasons, each file present in colvar folder is just a sub-sample of the original ones used in the paper.

***

# Run the sample
Once all the colvar files have been converted into binary ones, copy the `fit.x` program from `src` directory in the current one and run:
    `./fit.x < input`
