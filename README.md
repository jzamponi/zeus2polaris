# zeus2polaris
Convert from zeusTW dumpfiles into a POLARIS spherical grid format.

$ zeus2polaris --help
usage: zeus2polaris [-h] [-o OUTPUT] [-f {binary,ascii}]
                    [--constant-density CONSTANT_DENSITY] [--no-progress-bar]
                    [-v]
                    nframe

Convert from ZeusTW dumpfiles into a POLARIS binary grid

positional arguments:
  nframe                Frame number by which all input files are suffixed.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Name of the output grid. It will suffixed by the frame
                        number
  -f {binary,ascii}, --format {binary,ascii}
                        Format of the output grid. POLARIS uses binary but
                        ascii allows to check it.
  --constant-density CONSTANT_DENSITY
                        Option to create a constant density grid replacing the
                        input structure by a given number.
  --no-progress-bar     Disables the printing of a progress bar to standard
                        output in verbose mode.
  -v, --verbose
