This folder is for producing memory comparison as shown in Fig 4b in Section V.A.

How to run:
- Edit file `run_it.sh`:
    + Change `methodName` to either `MERIT` or `FastGradient`
    + (optionally) Set list of `N`
- Then run
```
cd <to_this_directory>
./run_it.sh
```
- Memory usage is shown in file `logging_<N>_<methodName>.txt` at line (near the end of file) starting with
`
Maximum resident set size (kbytes): 
`

