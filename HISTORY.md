# Revision History for `symlib`

## Rev. 2.0.3 (GLWH)
- There were many small changes on my local copy that I committed, mostly changes in the
absolute tolerances. They were made so long ago I don't remember why they were made.

## Rev 2.0.2
- Fixed a bug in both Smith Normal Form routines where it would
  attempt to divide by zero when one of diagonals was set to zero
  during one of the row/column operations.
- Added an optional absolute tolerance to the minkowski reduce code.

## Rev 2.0.1
- At the last revision, there were some left over markers from merge conflicts. Removed those.
- added 'regular = "true" to parameter tag for atom_pos; it had been changed from 'pointer' to 'allocatable'
- removed 'atom_pos' from the unit test xml file
- All of the unit tests seem to be passing fine *except* for two routines that are used in uncle. One is an out of bounds error in 'find_site_equivalencies'. The other is a missing file of some sort: `ERROR: can't create representation for /Users/hart/codes/staging/symmetry.find_site_equivalencies.g/tests/find_site_equivalencies_nSites.2/nSites.out. File does not exist.`. Will have to sort out those tests later.

# 2.0.0 Major revision (Jan 2019)
- Changed some pointers to be `allocatable arrays` instead. This was to keep `symlib` compatible with `enumlib` which used allocatable arrays and general avoids pointers.
- This counts as a major revision because it is a *breaking* change---the calling interface for several routines that are called from other programs are different now and the calling programs will have to change as a result.

## Revision 1.1.4
- Added a long int version of the `SmithNormalForm` algorithm.

## Revision 1.1.3
- Rewrote the `SmithNormalForm` algorithm to use a different pivoting
  strategy.
- Added an optional output err_ flag to the `SmithNormalForm`
  alforithm to report an error but not stop if the code overflows.

## Revision 1.1.2
- Changed the `bing_into_cell` subroutine to use `MOD` instead of a `do while` loop so improve runtimes.

## Revision 1.1.1
- Fixed bug (reported in [#19](https://github.com/msg-byu/symlib/issues/19)) where `max` and not `maxval` was being used when comparing arrays.

## Revision 1.1.0

- Added classes.f90 and intertolls.f90 to symlib from the polya repo
  (https://github.com/rosenbrockc/polya). Added them to the libsym.a
  compiled object.

## Revision 1.0.0

In rational_mathematics.f90, changed the comment for the SNF routine. The comments
on the calling interface were incorrect (A and M had been backwards).

In symmetry.f90, added some "regular = true" flags in the xml comments.

## Revision 0.0.1

Added unit tests for combinatorics.f90, it now has working unit tests for 66% of it's subroutines,
and its supporting combinatorics.xml file. Also cleaned up the combinatorics.f90 file and adde XML
documentation.

Added unit tests for compare_structures.f90, it has working unit tests for 100% of the subroutines
used in the UNCLE code. There are 5 other subroutines that are not used that have not yet been unit
tested. Also cleaned up and added XML documentation to the compare_structures.f90 file. Added the
compare_structures.xml file as well.

Added unit tests for `numerical_utilities.f90`, it has working unit tests for 90% of its
subroutines. Also changed the `equal` subroutines so that if the values are the same, it returns `true` no
matter the tolerance that is used. Added the XML documentation to `numerical_utilities.f90`, the unit
test files, and the numerical_utililties.xml file.

Added unit tests for `rational_mathematics.f90`, it has working unit tests for 60% of its
subroutines. The last 37% are waiting on a `fortpy` bug to get sorted out.

Added unit tests for `symmetry.f90`, but it is still a bit messy. For the 10, subroutines only 2 have fully functional unit tests. The rest are rank3 arrays, or require logical inputs. There are a few, marked in the xml doc, that fail because the saved and new files have different precision apparently `dp8` and `dp4`.

Added unit tests for `utilities.f90`, none work because of fortpy issues, but 100% of the subroutines have tests written.

Added unit tests for `vector_utilities.f90`. It is now 70% unit tested though the unit tests have some buggy floating-point issues.

## Initial Repository: Revision 0.0.0

The first few commits were to get the repo up to scratch and nice and clean with installation
instructions, etc. This includes the commit for new revision number 0.0.0. It includes an update of
the `*.xml` files defining the unit tests so that they work with the distribution directory's
structure.
