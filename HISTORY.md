# Revision History for `symlib`

## Revision 0.0.1

Added unit tests for combinatorics.f90, it now has working unit tests for 66% of it's subroutines, and it's supporting combinatorics.xml file. Also cleaned up the combinatorics.f90 file and adde XML documentation.

Added unit tests for compare_structures.f90, it has working unit tests for 100% of the subroutines used in the UNCLE code. There are 5 other subroutines that are not used that have not yet been unit tested. Also cleaned up and added XML documentation to the compare_structures.f90 file. Added the compare_structures.xml file as well.

Added unit tests for numerical_utilities.f90, it has working unit tests for 90% of it's subroutines. Also changed the equal subroutines so that if the value are the same it returns true no matter the tolerance that is used. Added the XML documentation to numerical_utilities.f90, the unit test files, and the numerical_utililties.xml file.

## Initial Repository: Revision 0.0.0

The first few commits were to get the repo up to scratch and nice and clean with installation instructions etc. This includes the commit for new revision number 0.0.0. It includes an update of the `*.xml` files defining the unit tests so that they work with the distribution directory's structure.