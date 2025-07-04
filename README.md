# T2Well_post

A script for post processing T2 well output simulations.
The script reads and parses the FStatus, FFlow, COFT and FOFT files and produces plots of the simulation results.


  ## Installation

  ### WINDOWS:
  1. **Install GIT**  
    Ensure GIT is installed on your machine. [GIT download link](https://accessit.equinor.com/Search/Search?term=%22GIT%22)

  2. **Install uv**  
    Open PowerShell or Command Prompt and run:  
    ```shell
    winget install --id=astral-sh.uv -e
    ```
  
  ### LINUX
  2. **Install uv**  
    Open PowerShell or Command Prompt and run:  
    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```
  
  3. **Install Python 3.11**  
    ```shell
    uv python install 3.11
    ```

  4. **Set up your project directory**  
    - Navigate to `C:\Appl` (or `C:\Appl\dev`):  
      ```shell
      cd \Appl
      ```
    - If the folder does not exist, create it:  
      ```shell
      mkdir Appl
      cd Appl
      mkdir dev  # optional
      ```

  5. **Create a new project directory**  
    ```shell
    mkdir myproject
    cd myproject
    ```

  6. **Initialize the environment**  
    ```shell
    uv init --python 3.11
    ```

  7. **Install this repository as a local package**  
    ```shell
    uv add git+https://github.com/equinor/T2Well_post
    ```

  8. **Upgrade to the latest version (if needed)**  
    ```shell
    uv add git+https://github.com/equinor/T2Well_post --upgrade
    ```


## Usage

Make sure that all tge input and output files are located in the same directory. this involves both the resilt reporting files (Fstatus, FFlow, COFT) and the system files (e.g. VERS).
By typing the line below, the script will produce 4 figures (1 per file), including all variables reported in the files. For the specific case of COFT and FOFT, the figures will include all items (mesh elements and connections).

```shell
#short form
uv run plot_T2output.py -p /path/to/input/file

#long form
uv run plot_T2output.py --input_file /path/to/input/file
```

If the user wants to plot fewer variables or fewer items you can type commands such as the one below:

```shell
#short form
uv run plot_T2output.py -p /path/to/input/file -fst -fst_var 1 3 4 -ffl -c -c_var 1 3 -f -f_item 4 5

#long form
uv run plot_T2output.py --input_file /path/to/input/file --FStatus -FStatus_vectors 1 3 4 --FFlow --COFT --COFT_vectors 1 3 --FOFT -FOFT_elements 4 5
```

For this particular example, note the inclusion of the letter “i" after the FOFT line. Such letter can be used for both COFT and FOFT. Anything before it refers to the variables included in the file. Anything after refers to the elements or connections.

In that line, the script will query the item/variables and plot as follow:

- Fstatus: 1st, 3rd and 4th variables
- FFlow: all variables
- COFT: 1st and 3rd variables for all connections reported in the COFT section
- FOFT: All variables for the 4th and 5th  mesh elements queried in the FOFT section

Maybe it is getting to complex and we can think of better ways of calling the script. Maybe including and input file.

| Command                             | Description                                             |
| ----------------------------------- | ------------------------------------------------------- |
| `-h`, `--help`                  | show help message and exit                              |
| `-logx`, `--log_scale_x`        | declare if a logarithmic X scale is required            |
| `-logy`, `--log_scale_y`        | declare if a logarithmic Y scale is required (COFT and FOFT only)            |
| `-fst`, `--FStatus`             | declare if FStatus is included                          |
| `-fst_var`, `--FStatus_vectors` | list which vectors of FStatus are included              |
| `-ffl`, `--FFlow`               | declare if FFlow is included                            |
| `-ffl_var`,` --FFlow_vectors`   | list which vectors of FFlow are included                |
| `-c`, `--COFT`                  | declare if COFT is included                             |
| `-c_var`, `--COFT_vectors`      | list which vectors of COFT are included                 |
| `-c_item`, `--COFT_connections` | list which connections of COFT are included             |
| `-f`, `--FOFT`                  | declare if FOFT is included                             |
| `-f_var` , `--FOFT_vectors`     | list which vectors of FOFT are included                 |
| `-f_item`, `--FOFT_elements`    | list which connections of FOFT are included             |
| `-xls`, `--print_xls`           | declare if output shall be printed in spreadsheet       |
| `-ts`, `--time_steps`           | generate plot of time steps                             |
| `-pcm`, `--pcolormesh`          | option to create a pcolormesh in FFlow and Fstatus plot |
| `-f_PT`, `--FOFT_PT`            | option to create a PT plot of FOFT elements             |
