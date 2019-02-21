# ACTIN 1.3
### Activity Indices Calculator

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00667/status.svg)](https://doi.org/10.21105/joss.00667)

Reads fits files from HARPS, HARPS-N and ESPRESSO spectrographs, rdb tables, and outputs user defined spectral activity indices (along with other relevant data).

Beta

### Requires the following Python modules:
- numpy
- matplotlib
- astropy
- appdirs (included)


### Installation:

Clone the github repository (in releases) to a directory of your choice and install via `python setup.py install`.


### Configuration file:

The `config_lines.txt` file is the line configuration file (instructions inside). This file is used to add line parameters to calculate any index as long as the line cores and bandpasses are inside the spectral range and spectral orders range (for 2d spectra) of the spectrograph. ACTIN will check this at start and give an error message if line parameters don't match the spectra.

This file is available from the directory each OS uses for storing user data (*)

To get your path to the configuration file call `actin` without any arguments.

The configuration file can be copied to another directory, modified, and used via `actin -cf dir/filename`.

NOTE: If not installed via pip, use `python actin.py` instead of `actin`.


### Quick start:

Usage:

`actin -h [help] -f [files_list] -i [indices_list] -rv [rv_list] -cf [config_file] -s [output_path] -lp [output_path/same/show] -obj [object_name] -tl [target_list] -del [True/False] -t [test_file_type] -frc [True/False]`


Arguments:

`-h` : Gives a description of the arguments available.

`-f` : List of files (formats S1D, S2D, e2ds, s1d, s1d_*_rv, or ADP) or rdb table(s) with headers `obj`, `date`, `bjd`, `wave`, `flux`, `error_pixel` (optional) to be read.

`-i` : List of indices to calculate. Indices ids must match the ones in the config file `config_lines.txt`.

`-rv` : List of RV values to calibrate wavelength. If not used, RVs are used from CCF files if available.

`-cf` : Path to configuration file. If not given the configuration file is read from installation directory To know installation directory call `actin` without arguments.

`-s` : Save output to .rdb table in specified path.

`-lp` : Save plots of the lines used to calculate the indices in the specified path. If `same` uses the path as specified in `-s` If `show` shows the plots without saving, useful to analyse the lines in one spectra.

`-tl` : List of stars to select from `files`.

`-del` : If `True` deletes any output file (data and logs only; only files that match current ACTIN call will be deleted) before reading the file list and saving output.

`-obj` : Object name to override the one from fits files in case the star has multiple names in different files (ex. Proxima, ProximaCen, Gl551). DO NOT USE WHEN READING FILES FROM MULTIPLE STARS.

`-t` : Tests the program using the test files. Options are `S2D`, `S1D`, `e2ds`, `s1d`, `adp` or `rdb` to test these type of files.

`-frc` : Use fractional pixels if `True` (default), use integral pixels if `False`. If using `False` and calculating the I_CaII index as given in the original config_lines.txt file ACTIN will simulate the values of 's_raw' from the HARPS pipeline. Note however that this option might induce artificial variations in the indices due to the use of integral pixels in the bandpasses.

#### Important:

When running ACTIN for a second time with the same data on the same output directory use `-del True` otherwise the program will detect the same dates,  ignore the measurements and not give any output.

When arguments accept lists, they can be given in the command line, e.g. `-tl Gl273 Gl581`, or from an ASCII file by using, e.g. `-tl $(cat target_list.txt)` where `target_list.txt` is a file with one column with the rows `Gl273` and `Gl581`.

### Testing the code with minimum arguments:

The example below will test the code using the test files provided in the package.

`actin -t e2ds`


### Example for multiple files:

`actin -f ../fits/*/*e2ds_A.fits -i I_CaII I_Ha -s ../output -del True -tl Gl273 Gl581`

This will execute ACTIN for all the subdirectories inside `../fits/` with files ending with `e2ds_A.fits`, calculate the indices `I_CaII` and `I_Ha`, output the data to `../output/star_names`, and, before running the code, delete any output file that was previously there, in this case `Gl273_HARPS_e2ds_data.rdb` and `Gl581_HARPS_e2ds_actin.rdb`. Only fits files belonging to the stars chosen in `-tl` will be read, in this case `Gl273` and `Gl581`. Since `-frc` is True by default, fractional pixels will be used to compute the indices.

---

(*) For OSX: `~/Library/Application Support/<AppName>`

For Windows: `C:\Documents and Settings\<User>\Application Data\Local Settings\<AppAuthor>\<AppName>` or possibly `C:\Documents and Settings\<User>\Application Data\<AppAuthor>\<AppName>`

For Linux: `~/.local/share/<AppName>`

---

For issues, problems, support or if you would like to contribute to the development of this code please contact Joao.Silva(at)astro.up.pt
