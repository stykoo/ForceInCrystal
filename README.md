# ForceInCrystal

ForceInCrystal is a software to investigate the behavior of a 2d
hexagonal lattice of Brownian particles when one of them is submitted
to an external force.


## Dependancies
* [Boost](http://www.boost.org/) for program options and unit testing.
* [SFML](http://www.sfml-dev.org/) for on-the-fly visualization of particles
in 2d

## Compilation
For compiling the source code, your compiler should support
[C++14](https://en.wikipedia.org/wiki/C%2B%2B14) (both
[g++](https://gcc.gnu.org/) and [clang++](http://clang.llvm.org/) do).

The compilation was only tested on Linux. It should in theory also work
on Windows and MacOS.

## Use of the program and documentation
The command-line arguments to the program are explained when
using the `--help` argument.

The code is documented using
[Doxygen](http://www.stack.nl/~dimitri/doxygen/)-compatible comments.
One can generate the html and latex documentation using `doxygen Doxyfile`.

## License
This program is free software and is released under
[GNU GPL](https://www.gnu.org/licenses/quick-guide-gplv3.html).
Feel free to contribute to the code!

* Boost is released under
[Boost Software License](http://www.boost.org/users/license.html).
* SFML is released by Laurent Gomila under
[zlib/png licence](http://www.sfml-dev.org/license.php).
