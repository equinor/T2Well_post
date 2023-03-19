# T2Well_post
A script for post processing T2 well output simulations.
The script reads and parses the FStatus, FFlow, COFT and FOFT files and produces plots of the simulation results.

## Installation:

Download/clone the contents of this repository

- Create a virtual environment
  ```shell
  python -m venv t2_venv
  ```
- Activate it
  ```shell
  . t2_venv/bin/activate
  ```
- Install requirements at location where the repository files are (unzipped).
  ```shell
  python -m pip install -r requirements.txt
  ```

## Running the script

By typing the line below, the script will produce 4 figures (1 per file), including all variables reported in the files. For the specific case of COFT and FOFT, the figures will include all items (mesh elements and connections).

```shell
python plot_T2output.py /path/to/input/file
```

If the user wants to plot fewer variables or fewer items you can type commands such as the one below:
```shell
python plot_T2output.py /path/to/input/file FStatus,1,3,4 FFlow COFT,1,3 FOFT,i,4,5
```

For this particular example, note the inclusion of the letter “i" after the FOFT line. Such letter can be used for both COFT and FOFT. Anything before it refers to the variables included in the file. Anything after refers to the elements or connections.

In that line, the script will query the item/variables and plot as follow:
- Fstatus: 1st, 3rd and 4th variables
- FFlow: all variables
- COFT: 1st and 3rd variables for all connections reported in the COFT section
- FOFT: All variables for the 4th and 5th  mesh elements queried in the FOFT section

Maybe it is getting to complex and we can think of better ways of calling the script. Maybe including and input file.


