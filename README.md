# symlib
Symmetry-related routines for cluster expansion and other codes that rely on symmetries of lattices and crystals.


## Compiling `symlib`

`symlib` ships with a Makefile in the `src/` directory. Although it has several rules for different use cases, for a quickstart, just type:

```
cd src/
make F90=[gfortran|ifort]
```

This will generate libraries `libutils.a libsym.a libcomparestructs.a librational.a libcombinatorics.a`. For API references on these libraries (which have to still be generated using `fortpy`), check out the wiki page.

## Running the Unit Tests

`symlib` is ~10% unit tested. Tests can be run using [fortpy](https://github.com/rosenbrockc/fortpy). We recommend using a virtualenv to run all the unit tests.

```
pip install fortpy
cd symlib
runtests.py src/
```

This will compile a series of drivers to test the `symlib` codes in a directory at `staging/`. Any errors or warnings generated will be redirected to `stdout`.

## Troubleshooting

So far, no one has *ever* had any problems installing and using this code :D.
