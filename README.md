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

```shell
python plot_T2output.py /path/to/input/file
```



