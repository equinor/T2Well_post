# -*- coding: utf-8 -*-

import os
from os import fstat
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.cm import ScalarMappable

import sys
import pandas as pd
import numpy as np

import argparse


plt.style.use('seaborn-white')


rcParams['figure.dpi'] = 300
rcParams['image.cmap'] = 'rainbow'
rcParams['ytick.major.size'] = 3
rcParams['xtick.major.size'] = 3
rcParams['font.size'] = 6
rcParams['lines.linewidth'] = 0.75




#### Parsing functions

def read_FStatus(ip_file,eos):

    """
    Function to parse FStatus file.
    It creates a pandas dataframe
    """

    print('Processing {:s} file'.format(ip_file))

    with open(ip_file) as fs:
        next(fs)
        fs_labels = next(fs)
    
    fs_labels = fs_labels.split('=')[1]
    
    if eos == 'ECO2N':
        fs_labels = fs_labels.split()
        fs_header = fs_labels[2:]
    
    else:
        fs_labels = fs_labels.split(",")
        fs_labels = [item.split('(')[0] for item in fs_labels]
        fs_header = fs_labels[4:]


    fstatus = pd.read_csv(ip_file, skiprows=3, names=fs_labels)
    fstatus = fstatus.apply(pd.to_numeric, errors='coerce')
    fstatus = fstatus.fillna(0)
    fstatus['Pres']/= 1e5


    return fs_header, fstatus

def read_FFlow(ip_file, eos):

    """
    Function to parse FFlow file.
    It creates a pandas dataframe
    """


    print('Processing {:s} file'.format(ip_file))

    with open(ip_file) as ff:
        if eos != 'ECO2N':
            next(ff)
        ff_labels = next(ff)
    
    if eos == 'ECO2N':
        ff_labels = ff_labels.split()
        r_skip = 1

        ff_header = ff_labels[2:]
    
    else:
        ff_labels = ff_labels.split('=')[1]
        ff_labels = ff_labels.split(",")
        ff_labels = [item.split('(')[0] for item in ff_labels]
        ff_labels.append('V_mix')
        # add dummy header to column created at the end
        ff_labels.append('trail')
        r_skip = 3

        ff_header = ff_labels[4:]



    fflow = pd.read_csv(ip_file, skiprows=r_skip, names=ff_labels)
    fflow = fflow.apply(pd.to_numeric, errors='coerce')
    fflow = fflow.fillna(0)

    if eos != 'ECO2N':
        fflow = fflow.drop('trail', axis=1)
        ff_header.pop()


    return ff_header, fflow


def read_COFT(ip_file, eos):

    """
    Function to parse coft file.
    It creates a pandas dataframe
    """


    print('Processing {:s} file'.format(ip_file))


    coft = pd.read_csv(ip_file, header=None, index_col=0, low_memory=False)

    coft = coft.apply(pd.to_numeric, errors='coerce')

    #Drop empty columns
    coft = coft.dropna(axis=1, how='all')
    
    #Replace remaining NaN values for zero
    coft = coft.fillna(0)
    
    if eos == 'ECO2N':
        coft_var = ['FLO(GAS)', 'FLO(aq.)', 'VEL(GAS)', 'VEL(LIQ.)', 'FLO(NaCl)' , 'FLO(CO2)' ,'FLOH']
        c_idx_df = coft.loc[:,2::8]
    else:
        coft_var = ['FLO(aq.)', 'FLO(LIQ.)', 'FLO(GAS)', 'FLO(CO2)']
        c_idx_df = coft.loc[:,2::5]



    #Retrieve connection indexes and delete columns
    

    c_idx = c_idx_df.drop_duplicates().values.flatten()
    coft = coft.drop(columns=c_idx_df.columns)


    print('{:4d} connections reported in coft.'.format(len(c_idx)))


    #Define multiIndex to rename columns
    coft_col = pd.MultiIndex.from_product([c_idx,coft_var])
    coft_col = coft_col.insert(loc=0, item='time')


    coft.columns = coft_col


    return coft_var, c_idx, coft

def read_FOFT(ip_file, eos):

    """
    Function to parse foft file.
    It creates a pandas dataframe
    """


    print('Processing {:s} file'.format(ip_file))


    foft = pd.read_csv(ip_file, header=None, index_col=0, low_memory=False)
    
    foft = foft.apply(pd.to_numeric, errors='coerce')
    
    #Drop empty columns
    foft = foft.dropna(axis=1, how='all')
    
    #Replace remaining NaN values for zero
    foft = foft.fillna(0)
    
    #Define headers
    if eos == 'ECO2N':
        foft_var = ['Pres', 'Sg', 'XNACL', 'XCO2liq', 'T']
    else:
        foft_var = ['Pres', 'T', 'S_liquid', 'S_gas', 'XCO2liq']
    
    #Retrieve element indexes and delete columns
    e_idx_df = foft.loc[:,2::6]


    e_idx = e_idx_df.drop_duplicates().values.flatten()

    
    print('{:4d} grid elements reported in foft.'.format(len(e_idx)))

    foft = foft.drop(columns=e_idx_df.columns)


    #Define multiIndex to rename columns
    foft_col = pd.MultiIndex.from_product([e_idx,foft_var])
    foft_col = foft_col.insert(loc=0, item='time')


    foft.columns = foft_col
    
    foft.iloc[:,1::5] = foft.iloc[:,1::5]/1e5

    



    return foft_var, e_idx, foft

def read_ipMESH(fname):

    """
    Parse TOUGH input file and extract grid parameters
    from ELEME and CONNE sections into a dataframe
    """

    mesh = dict()

        
    with open(fname, 'r') as ip_file:
        f = ip_file.read()
        
        if 'ELEME' in f:
            read_MESH_file = False
        else:
            read_MESH_file = True


    if read_MESH_file:
        print("Read ELEME and CONNE from MESH file")
        with open(r'MESH') as mesh_file:
            f = mesh_file.read()

    ELEME_i = f.find('ELEME')
    ELEME_f = ELEME_i + f[ELEME_i:].find('CONNE')
    ELEME_raw = f[ELEME_i:ELEME_f].split('\n')[1:]


    ELEME = pd.DataFrame(data = ELEME_raw, index = 1+np.arange(len(ELEME_raw)))
    ELEME = ELEME.drop(index=ELEME[ELEME[0].str.strip().str.len()==0].index)
    

    ELEME['ElName'] = ELEME[0].str[:5]
    ELEME['NSEQ'] = ELEME[0].str[5:10]
    ELEME['NADD'] = ELEME[0].str[10:15]
    ELEME['MAT'] = ELEME[0].str[15:20]
    ELEME['VOLX'] = pd.to_numeric(ELEME[0].str[20:30], 'coerce')
    ELEME['AHTX'] = pd.to_numeric(ELEME[0].str[30:40], 'coerce')
    ELEME['PMX'] = pd.to_numeric(ELEME[0].str[40:50], 'coerce')
    ELEME['X'] = pd.to_numeric(ELEME[0].str[50:60], 'coerce')
    ELEME['Y'] = pd.to_numeric(ELEME[0].str[60:70], 'coerce')
    ELEME['Z'] = pd.to_numeric(ELEME[0].str[70:80], 'coerce')

    #Remove zero volume elements set for Dirichlet condition
    ELEME.loc[ELEME.VOLX.isna(),'VOLX'] = 0

    #Reset index to preserve correct index values
    ELEME  = ELEME.loc[~(ELEME.VOLX==0)].copy().reset_index(drop=True)
    ELEME.index = ELEME.index+1
    
    ELEME = ELEME.drop(columns=0)


    CONNE_i = f.find('CONNE')
    
    if read_MESH_file:
        CONNE_f = CONNE_i + f[CONNE_i:].find('+++')
        CONNE_raw = f[CONNE_i:CONNE_f].split('\n')[1:-1]
    else:
        f = f.replace('\n \n', '\n\n')
        CONNE_f = CONNE_i + f[CONNE_i:].find('\n\n')
        CONNE_raw = f[CONNE_i:CONNE_f].split('\n')[1:]

    CONNE = pd.DataFrame(data = CONNE_raw, index = 1+np.arange(len(CONNE_raw)))
    CONNE['EL1'] = CONNE[0].str[:5]
    CONNE['EL2'] = CONNE[0].str[5:10]
    CONNE['NSEQ'] = CONNE[0].str[10:15]
    CONNE['NAD1'] = CONNE[0].str[15:20]
    CONNE['NAD2'] = CONNE[0].str[20:25]
    CONNE['ISOT'] = pd.to_numeric(CONNE[0].str[25:30], 'coerce')
    CONNE['D1'] = pd.to_numeric(CONNE[0].str[30:40], 'coerce')
    CONNE['D2'] = pd.to_numeric(CONNE[0].str[40:50], 'coerce')
    CONNE['AREAX'] = pd.to_numeric(CONNE[0].str[50:60], 'coerce')
    CONNE['BETAX'] = pd.to_numeric(CONNE[0].str[60:70], 'coerce')
    CONNE['SIGX'] = pd.to_numeric(CONNE[0].str[70:80], 'coerce')

    CONNE = CONNE.drop(columns=0)


    return ELEME, CONNE




#Plotting functions

name_variable = ['Depth', 'Dis', 'cumDepth', 'Dgas', 'D_aqueous' , 'D_liquid', 'D_gas', 'FLO(CO2)', 'FLO(GAS)', 
                'FLO(NaCl)', 'FLO(aq.)', 'FLO(LIQ.)',  'q_aqueous', 'q_liquid', 'q_gas' , 'FLOH', 'Fgas', 
                'Fliq', 'Pres', 'Sg', 'S_aqueous', 'S_liquid', 'S_gas', 'T', 'Time', 
                'Umix', 'VEL(GAS)', 'VEL(LIQ.)', 'VGas',  'V_aqueous', 'V_liquid', 'V_gas', 'V_mix',
                'VLiq', 'XCO2liq', 'XNACL', 'var1', 'var2', 'var3', 'var4', 'WellID']

unit_variable = ['m']*3 +    4*['kg/m$^3$'] + 8*['kg/s'] + ['W'] + 2*['kg/s'] + ['bar'] + ['m$^3$/m$^3$']*4 + ['$\degree$C', 's'] + 9*['m/s'] + 2*['kg/kg'] + 5*['-']
unit_variable_v2 = ['m']*3 + 4*['kg/m3']    + 8*['kg/s'] + ['W'] + 2*['kg/s'] + ['bar'] + ['m3/m3'      ]*4 + ['degC', 's']       + 9*['m/s'] + 2*['kg/kg'] + 5*['-']

unit_lims = {'Sg':(0,1)}

units_dict = dict(zip(name_variable,unit_variable))
units_dict_v2 = dict(zip(name_variable,unit_variable_v2))

perm_dict = {-1:'-x', 1:'+x', -2:'-y', 2:'+y', -3:'-z', 3:'+z'}


def secondary_scale(log_bool, ax):
    """
    Sets a secondary scale with human readable units
    Takes two parameters
    log_bool: Boolean set as True if the scale is logarithmic, False if linear
    ax: The matploltib artist axis

    Returns the input axis with a new secosndary x-axis on top
    """
    if log_bool:
        ax.set_xscale('log')


        sec_x = ax.twiny()

        sec_x.set_xscale('log')
        x_ticks = np.array([1,60,3600,3600*24,3600*24*30, 3600*24*365.25, 3600*24*365.25*100, 3600*24*365.25*1000])
        x_labels = np.array(['1 s', '1 min', '1 h', '1 day', '1 month', '1 yr', '100 yr', '1 kyr'])
        sec_x.set_xticks(x_ticks)
        sec_x.set_xticklabels(x_labels)
        sec_x.set_xlim(ax.get_xlim())
        sec_x.set_xlabel('time')
            
    else:
        
        
        t_min, t_max = ax.get_xlim()
        
        scales = [60, 3600, 3600*24, 3600*24*365.25, 3600*24*365.25*1000]
        units = ['min', 'hr', 'd', 'yr', 'kyr']



        n_ticks = 5


        fix = np.arange(n_ticks+1)
            
        sec_x = ax.twiny()
        for s,u in zip(scales, units):
            scaled_tmax = t_max/s
            
            order = np.log10(scaled_tmax)


            if order>0 and order<=3:
                if order<2:
                        if scaled_tmax/n_ticks<1:

                                step = (5*(1+10*(scaled_tmax/n_ticks))//5)/10
                                x_labels = fix*step



                        elif scaled_tmax<=10:
                                x_labels = np.arange(np.ceil(scaled_tmax))
                        else:
                                step = scaled_tmax//n_ticks
                                x_labels = np.arange(np.ceil(scaled_tmax), step = step)
                
                else:
                        step = 50*(1+(scaled_tmax//n_ticks)//50)
                        x_labels = fix*step
                
                x_labels = np.round(x_labels, 1)
                x_ticks = x_labels*s
                sec_unit = u



                                        
                sec_x.set_xticks(x_ticks)
                sec_x.set_xticklabels(x_labels)
                sec_x.set_xlim(ax.get_xlim())
                sec_x.set_xlabel('time [{:s}]'.format(sec_unit))    



def plot_Ffigure(title,df,df_vars, logscale, EOS):
    """
    Function for plotting the FFlow or FStatus files
    Takes three parameters
    title: FStatus or FFlow
    df: the dataframe corresponding to each file
    vars: the variables set for plotting
    logscale: boolean if the plot is needed in logscale
    EOS: Specifics about the naming of variables

    Returns and save a figure
    """
    #filter index bars
    index_var = ['WellID',  'Time',  'Dis',  'cumDepth', 'Depth']
    plot_vars = [item for item in df_vars if item not in index_var]
    
    #Set plot dimensions
    w = 5
    h = (2.5-1)*len(plot_vars)
    rcParams['figure.figsize'] = [w,h]

    #Create axes
    fig, axs = plt.subplots(len(plot_vars),1, sharex=True)


    #Pick right depth vector depending on labels difference in EOS version
    if EOS == 'ECO2N':
        Z_label = 'Depth'
    else:
        Z_label = 'Dis'
    
    #Define x, y dimensions for creating contour plot
    size_j = df[Z_label].drop_duplicates().shape[0]
    size_i = int(df['Time'].shape[0]/size_j)

    #Define the X and Y axis
    X = df['Time'].values.reshape(size_i, size_j)
    Y = df[Z_label].values.reshape(size_i, size_j)


    for var_idx, var in enumerate(plot_vars):

        #Define limits for saturation variables
        if var in ['Sg', 'S_aqueous', 'S_liquid', 'S_gas']:
            unit_lims[var] = (0,1)


        if len(plot_vars)>1:
            ax = axs[var_idx]
        else:
            ax = axs


        Z = df[var].values.reshape(size_i, size_j)
        

        try:
            cf = ax.contourf(X, Y, Z, vmin = unit_lims[var][0], vmax = unit_lims[var][1])
            cb = plt.colorbar(ScalarMappable(norm=cf.norm, cmap=cf.cmap), ax=ax,
                                ticks=np.linspace(unit_lims[var][0], unit_lims[var][1], 6))
        except:
            cf = ax.contourf(X, Y, Z)
            cb = plt.colorbar(cf, ax=ax)
        
        cb.set_label('{:s} [{:s}]'.format(var,units_dict[var]))
        ax.set_ylabel('depth [m]')

        ax.invert_yaxis()

        if var == plot_vars[-1]:
            ax.set_xlabel('time [s]')

        if var == plot_vars[0]:
            secondary_scale(logscale, ax)

    fig.suptitle(title)
    fig.align_ylabels()
    fig.tight_layout(rect=[0,0,1,0.98])
    fig.savefig(r'fig_{:s}.png'.format(title))




def plot_OFT(title, df, items, df_vars, logscale, mesh_eleme, mesh_conne):
    """
    Function for plotting the COFT or FOFT files
    Takes four parameters
    title: COFT or FOFT
    df: the dataframe corresponding to each file
    items: The connection(s) or element(s) set for plottinh
    df_vars: the variables set for plotting

    logscale: boolean if the plot is needed in logscale
    mesh_eleme: table with grid elements
    mesl_conne: table with grid connections

    Returns and save a figure
    """
    w = 6
    h = (2.8-1)*len(df_vars)
    
    w = 8.1
    h = 5.85
    rcParams['figure.figsize'] = [w,h]

    mesh_eleme_v2 = mesh_eleme.copy()
    mesh_eleme_v2 = mesh_eleme_v2.set_index(('ElName'))


    fig, axs = plt.subplots(len(df_vars),1, sharex=True)

    for var_idx, var in enumerate(df_vars):

        print(var_idx, var)


        if len(df_vars)>1:
            ax = axs[var_idx]
        else:
            ax = axs

        for item in items:
            # print('plot item {:d} in {:s} plot'.format(item,var))

            if title =='coft':
   
                el1 = mesh_conne.loc[item, 'EL1']
                el2 = mesh_conne.loc[item, 'EL2']
                k_dir = mesh_conne.loc[item, 'ISOT']

                mat1 = mesh_eleme_v2.loc[el1,'MAT']
                mat2 = mesh_eleme_v2.loc[el2,'MAT']
                item_label =  '{:s}>{:s}({:s} to {:s} in {:s} dir.)'.format(el1,el2,mat1, mat2, perm_dict[k_dir])



            elif title =='foft':
                el = mesh_eleme.loc[item,'ElName']
                mat = mesh_eleme.loc[item,'MAT']
                item_label = '{:<4d}{:s}({:s})'.format(item,el,mat)
            df.plot(x='time', y=(item, var), ax=ax, label=item_label, legend = False)

        if var == df_vars[-1]:
            ax.set_xlabel('time [s]')

        if var == df_vars[0]:
            secondary_scale(logscale, ax)
        

        
        
        ax.set_ylabel('{:s} [{:s}]'.format(var, units_dict[var]))
        
    if title == 'foft':
        ax.set_yscale('log')

    fig.suptitle(title)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right')
    fig.align_ylabels()
    fig.tight_layout(rect=[0,0,0.75,0.98])
    fig.savefig(r'fig_{:s}.png'.format(title))






#Define variables for plotting based on input arguments

def get_EOS(fname):
    """Retrieves EOS from output file"""
    with open(fname, 'r', encoding='latin1') as of:
        lines = of.readlines()

    for line in lines:
        line = line.strip()
        if line.startswith('EOS'):
            EOS = line.split()[5].strip('*')
            return EOS

def plot_specs(ip_args, plot_bool, files):
    """Defines the files and variables that will be plotted"""   
    plot_dict = dict()

    if len(ip_args)>2:
        
        #Turn off all plotting
        for f in plot_bool:
            plot_bool[f] = False


        #Based on arguments turn on the selected files to plot

        for arg in ip_args[2:]:
        
            arg = arg.split('i')
            arg_v = arg[0].strip(",").split(",")

            try:
                arg_i = arg[1].strip(",").split(",")
            except:
                arg_i = None

            print('file queried:', arg_v)
            print('items queried:', arg_i)

            #Check if FStatus is being queried
            if any(item.lower().startswith('fstatus') for item in files) and arg_v[0].lower() == r'fstatus':
                plot_bool[r'fstatus'] = True
                if len(arg_v)>1:
                    plot_dict[r'fstatus'] = arg_v[1:]
                else:
                    plot_dict[r'fstatus'] = 'all'
                    

            #Check if FFLow is being queried
            if any(item.lower().startswith('fflow') for item in files) and arg_v[0].lower() == r'fflow':

                plot_bool[r'fflow'] = True

                if len(arg_v)>1:
                    plot_dict[r'fflow'] = arg_v[1:]
                else:
                    plot_dict[r'fflow'] = 'all'


            #Check if coft is being queried
            elif any(item.startswith('COFT') for item in files) and arg_v[0].lower() == r'coft':
                plot_bool[r'coft'] = True
                plot_dict[r'coft'] = dict()
                if len(arg_v)>1:
                    plot_dict[r'coft']['var'] = arg_v[1:]
                else:
                    plot_dict[r'coft']['var'] = 'all'

                if arg_i is None:
                    plot_dict[r'coft']['item'] = 'all'
                else:
                    plot_dict[r'coft']['item'] = arg_i
                

            #Check if foft is being queried
            elif any(item.startswith('FOFT') for item in files) and arg_v[0].lower() == r'foft':

                plot_bool[r'foft'] = True
                plot_dict[r'foft'] = dict()
                if len(arg_v)>1:
                    plot_dict[r'foft']['var'] = arg_v[1:]
                else:
                    plot_dict[r'foft']['var'] = 'all'

                if arg_i is None:
                    plot_dict[r'foft']['item'] = 'all'
                else:
                    plot_dict[r'foft']['item'] = arg_i

        


    else:

        for f in plot_bool:

            plot_dict[f] = 'all'

            if f in ['coft', 'foft']:
                plot_dict[f]=dict()
                plot_dict[f]['var'] = 'all'
                plot_dict[f]['item'] = 'all'
    

    return plot_bool, plot_dict


def Excel_printer(fnames_map, ip_file, eleme, conne):
    """Orchestrates the printing of tables on a spreadsheet"""
    spreadsheet = ip_file.split(".")[0]+".xlsx"
    pd_units = pd.Series(units_dict_v2)
    
    if r'fstatus' in fnames_map.keys():
        #Add column names to FStatus

        fs_var, fs = read_FStatus(fnames_map['fstatus'], EOS)
        fs_row1 = fs.columns.to_list()
        fs_row2 = pd_units[fs_row1].to_list()

        fs_cols = pd.MultiIndex.from_tuples(list(zip( fs_row1, fs_row2)))

        fs.columns = fs_cols


    if r'fflow' in fnames_map.keys():

        #Add column names to FFlow

        ff_var, ff = read_FFlow(fnames_map['fflow'], EOS)
        ff_row1 = ff.columns.to_list()
        ff_row2 = pd_units[ff_row1].to_list()

        ff_cols = pd.MultiIndex.from_tuples(list(zip( ff_row1, ff_row2)))

        ff.columns = ff_cols


    if r'coft' in fnames_map.keys():

        #Add column names to coft
        coft_var, coft_idx, coft = read_COFT(fnames_map['coft'], EOS)

        eleme2 = eleme.copy()
        eleme2 = eleme2.set_index('ElName')


        query_x = conne['ISOT']==1
        query_y = conne['ISOT']==2
        query_z = conne['ISOT']==3





        coft_row1 = coft.columns.get_level_values(0).to_list()
        coft_row2 = (conne.loc[coft_row1[1:],'EL1']+conne.loc[coft_row1[1:],'EL2']).to_list()
        coft_row2 = [''] + coft_row2
        coft_row3 = coft.columns.get_level_values(1).to_list()
        coft_row4 = ['s'] + list(pd_units[coft.columns.get_level_values(1)[1:].to_list()].values)

        coft_cols = pd.MultiIndex.from_tuples(list(zip( coft_row1, 
                                                        coft_row2, 
                                                        coft_row3,
                                                        coft_row4)))

        coft.columns = coft_cols

    if r'foft' in fnames_map.keys():

        #Add column names to foft
        foft_var, foft_idx, foft = read_FOFT(fnames_map['foft'], EOS)

        foft_row1 = foft.columns.get_level_values(0).to_list()
        foft_row1[0] = 'cell_idx'

        foft_row2 = eleme.loc[foft_row1[1:],'ElName'].to_list()
        foft_row2 = ['cell_name'] + foft_row2

        foft_row2a = ['cell_X [m]'] + eleme.loc[foft_row1[1:], 'X'].to_list()
        foft_row2b = ['cell_Y [m]'] + eleme.loc[foft_row1[1:], 'Y'].to_list()
        foft_row2c = ['cell_Z [m]'] + eleme.loc[foft_row1[1:], 'Z'].to_list()
        foft_row2d = ['cell_mat'] + eleme.loc[foft_row1[1:], 'MAT'].to_list()

        foft_row3 = foft.columns.get_level_values(1).to_list()
        foft_row3[0] = 'time'
        foft_row4 = ['s'] + list(pd_units[foft.columns.get_level_values(1)[1:].to_list()].values)

        foft_cols = pd.MultiIndex.from_tuples(list(zip( foft_row1, 
                                                        foft_row2, 
                                                        foft_row2a, 
                                                        foft_row2b, 
                                                        foft_row2c, 
                                                        foft_row2d, 
                                                        foft_row3,
                                                        foft_row4)))

        foft.columns = foft_cols


    with pd.ExcelWriter(spreadsheet) as writer:
        if r'fstatus' in fnames_map.keys():
            fs.to_excel(writer, sheet_name=r'FStatus')
        if r'fflow' in fnames_map.keys():
            ff.to_excel(writer, sheet_name=r'FFlow')
        if r'coft' in fnames_map.keys():
            coft.to_excel(writer, sheet_name=r'coft')
        if r'foft' in fnames_map.keys():
            foft.to_excel(writer, sheet_name=r'foft')



def plotter_manager(fnames_map, plot_bool, eleme, conne):
    """Orchestrates the plotting of figures"""
    for ftype in fnames_map:    
        
        """debugging
        print(f'filetype {ftype} is named {fnames_map[ftype]}')
        if plot_bool[ftype]:
            print(f'{ftype} will be plotted')
        else:
            print(f'{ftype} will not be plotted')
        """
    
        file = fnames_map[ftype] #filename

        plot_f = plot_bool[ftype] #check if file will be included in plot


        if plot_f:
            print(f'Plotting {file} data\n')

            if ftype in ['fflow', 'fstatus']:
                queried_vars = plot_dict[ftype]

            

                df_vars, df = parse_dict[ftype](file, EOS)
                
                selected_var = df_vars


                if queried_vars != 'all':
                    vars_sel = []
                    for var in queried_vars:
                        try:
                            var_idx = (int(var)-1)%len(df_vars)
                            vars_sel.append(df_vars[var_idx])

                        except:
                            for var_name in df_vars:
                                if var.lower() == var_name.lower():
                                    vars_sel.append(var_name)
                    selected_var = vars_sel

                print(f'Plot for "{file}" includes: {" ".join(selected_var)}')

                plot_Ffigure(ftype, df,selected_var, logscale, EOS)




            elif ftype in ['coft', 'foft']:
                # print(plot_dict)

                queried_vars = plot_dict[ftype]['var']
                queried_items = plot_dict[ftype]['item']
                
                # print(queried_vars)
                
                
                df_vars, df_items, df = parse_dict[ftype](file, EOS)
                

                selected_var = df_vars
                selected_items = df_items

                # print(selected_items)

                if queried_vars != 'all':
                    vars_sel = []
                    for var in queried_vars:
                        try:
                            var_idx = (int(var)-1)%len(df_vars)
                            vars_sel.append(df_vars[var_idx])

                        except:
                            for var_name in df_vars:
                                if var.lower() == var_name.lower():
                                    vars_sel.append(var_name)
                    selected_var = vars_sel

                if queried_items != 'all':
                    item_sel = []
                    for item in queried_items:
                        try:
                            item_idx = (int(item)-1)%len(df_items)
                            item_sel.append(df_items[item_idx])

                        except:
                            for item_name in df_items:
                                if item.lower() == item_name.lower():
                                    item_sel.append(item_name)
                    selected_items = item_sel


                print('{:s} plotted variables include: {:s}'.format(file, ' '.join(selected_var)))
                print('{:s} plotted items include: {:s}'.format(file, ' '.join(map(str,selected_items))))
                plot_OFT(ftype, df, selected_items, selected_var, logscale, mesh_eleme=eleme, mesh_conne=conne)


if __name__ == '__main__':
    args = sys.argv
    print(args)
    logscale = False

    ip_path = args[1]
    ip_dirname = os.path.dirname(ip_path)

    old_path = os.getcwd()

    if len(ip_dirname)>0:
        os.chdir(ip_dirname)

    
    raw_names = tuple([r'fflow', r'fstatus', r'coft', r'foft'])
    fnames = []
    plot_bool = dict()
    fnames_map=dict()

    #Map out files linked to input file
    for f in os.listdir():
        if f.lower().startswith(raw_names):
            f_size = os.path.getsize(f)
            if f_size>0:
                # print(f'{f} file exists and will be added to fnames_map and plot_bool')
                fnames.append(f)
                flabel = f.lower().split('_')[0]
                fnames_map[flabel] = f
                plot_bool[flabel] = True

        elif f.endswith(tuple(['.in', '.inp'])):
            ip_file = f

        elif f.endswith('out'):
            op_file = f


    EOS = get_EOS(op_file)
    print(f'T2Well input file: {ip_file}')
    print(f'T2Well output file: {op_file}')
    print('EOS version: {:s}'.format(EOS))

    #Define horizontal scale type
    if 'log' in args:
        print('Plots will be plot in logarithmic scale')
        logscale = True
        args.remove('log')
    else:
        print('Plots will be plot in linear scale')


    #Define how files will be parsed
    parse_dict = dict()

    for file in fnames_map:
        if file.lower().startswith('fflow'):
            parse_dict[file]=read_FFlow
        elif file.lower().startswith('fstatus'):
            parse_dict[file]=read_FStatus
        elif file.startswith('coft'):
            parse_dict[file]=read_COFT
        elif file.startswith('foft'):
            parse_dict[file]=read_FOFT


    #Define which files and variables will be plotted
    plot_bool, plot_dict = plot_specs(args, plot_bool, fnames)
    

    print(f'plot_bool is {plot_bool}')
    print(f'plot_dict is {plot_dict}')

    eleme, conne = read_ipMESH(ip_file)

    
    plotter_manager(fnames_map, plot_bool, eleme, conne)

    #Query and index data

    # print(f'fnames_map is {fnames_map}')

    #Print data in spreadsheet
    print_Excel = True

    xls_output = input('Do you want to store output files as spreadsheet (Y/N). Default is Yes?\t')

    if xls_output.lower().startswith('n'):
        print_Excel = False


   
    if print_Excel:
        Excel_printer(fnames_map, ip_file, eleme, conne)