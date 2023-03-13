# -*- coding: utf-8 -*-

import os
from os import fstat
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.cm import ScalarMappable
import argparse


import sys
import pandas as pd
import numpy as np

plt.style.use('seaborn-white')


rcParams['figure.dpi'] = 300
rcParams['image.cmap'] = 'rainbow'
rcParams['ytick.major.size'] = 3
rcParams['xtick.major.size'] = 3
rcParams['font.size'] = 6
rcParams['lines.linewidth'] = 0.75






#### Parsing functions

def read_FStatus(ip_file):

    """
    Function to parse FStatus file.
    It creates a pandas dataframe
    """

    print('Processing {:s} file'.format(ip_file))

    with open(ip_file) as fs:
        next(fs)
        fs_header = next(fs)
    
    fs_header = fs_header.split('=')[1]
    fs_header = fs_header.split()

    fstatus = pd.read_csv(ip_file, skiprows=3, names=fs_header)
    fstatus = fstatus.apply(pd.to_numeric, errors='coerce')
    fstatus = fstatus.fillna(0)
    fstatus['Pres']/= 1e5


    return fs_header[2:], fstatus

def read_FFlow(ip_file):

    """
    Function to parse FFlow file.
    It creates a pandas dataframe
    """


    print('Processing {:s} file'.format(ip_file))

    with open(ip_file) as ff:
        ff_header = next(ff)
    
    ff_header = ff_header.split()


    fflow = pd.read_csv(ip_file, skiprows=1, names=ff_header)
    fflow = fflow.apply(pd.to_numeric, errors='coerce')
    fflow = fflow.fillna(0)


    return ff_header[2:], fflow

def read_COFT(ip_file):

    """
    Function to parse COFT file.
    It creates a pandas dataframe
    """


    print('Processing {:s} file'.format(ip_file))


    coft = pd.read_csv(ip_file, header=None, index_col=0, low_memory=False)

    coft = coft.apply(pd.to_numeric, errors='coerce')

    #Drop empty columns
    coft = coft.dropna(axis=1, how='all')
    
    #Replace remaining NaN values for zero
    coft = coft.fillna(0)
    
    
    coft_var = ['FLO(GAS)', 'FLO(aq.)', 'VEL(GAS)', 'VEL(LIQ.)', 'FLO(NaCl)' , 'FLO(CO2)' ,'FLOH']


    #Retrieve connection indexes and delete columns
    c_idx_df = coft.loc[:,2::8]

    c_idx = c_idx_df.drop_duplicates().values.flatten()
    coft = coft.drop(columns=c_idx_df.columns)


    print('{:4d} connections reported in COFT.'.format(len(c_idx)))


    #Define multiIndex to rename columns
    coft_col = pd.MultiIndex.from_product([c_idx,coft_var])
    coft_col = coft_col.insert(loc=0, item='time')


    coft.columns = coft_col


    return coft_var, c_idx, coft

def read_FOFT(ip_file):

    """
    Function to parse FOFT file.
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
    foft_var = ['Pres', 'Sg', 'XNACL', 'XCO2liq', 'T']
    
    #Retrieve element indexes and delete columns
    e_idx_df = foft.loc[:,2::6]


    e_idx = e_idx_df.drop_duplicates().values.flatten()

    
    print('{:4d} grid elements reported in FOFT.'.format(len(e_idx)))

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

    mesh['ELEME'] = ELEME
    mesh['CONNE'] = CONNE


    return mesh




#Plotting functions

name_variable = ['Depth', 'Dgas', 'FLO(CO2)', 'FLO(GAS)', 
                'FLO(NaCl)', 'FLO(aq.)', 'FLOH', 'Fgas', 
                'Fliq', 'Pres', 'Sg', 'T', 'Time', 
                'Umix', 'VEL(GAS)', 'VEL(LIQ.)', 'VGas', 
                'VLiq', 'XCO2liq', 'XNACL']

unit_variable = ['m', 'kg/m$^3$'] + 4*['kg/s'] + ['W'] + 2*['kg/s'] + ['bar', 'm$^3$/m$^3$', '$\degree$C', 's'] + 5*['m/s'] + 2*['kg/kg']

unit_variable_v2 = ['m', 'kg/m3'] + 4*['kg/s'] + ['W'] + 2*['kg/s'] + ['bar', 'm3/m3', 'degC', 's'] + 5*['m/s'] + 2*['kg/kg']

unit_lims = {'Sg':(0,1)}

units_dict = dict(zip(name_variable,unit_variable))
units_dict_v2 = dict(zip(name_variable,unit_variable_v2))

perm_dict = {-1:'-x', 1:'+x', -2:'-y', 2:'+y', -3:'-z', 3:'+z'}


def secondary_scale(log_bool, ax):
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



def plot_Ffigure(title,df,vars, logscale):
    
    w = 5
    h = (2.5-1)*len(vars)
    rcParams['figure.figsize'] = [w,h]

    fig, axs = plt.subplots(len(vars),1, sharex=True)

    size_j = df['Depth'].drop_duplicates().shape[0]
    size_i = int(df['Time'].shape[0]/size_j)

    df['Time_d'] = df['Time']/(3600*24)
    
    
    if title == r'FFlow':
        for var in vars[0:2]:
            unit_lims[var] = ( np.percentile(df[vars[0:2]], 0.25), np.percentile(df[vars[0:2]], 99.75) )
        for var in vars[2:]:
            unit_lims[var] = ( np.percentile(df[vars[2:]], 0.25), np.percentile(df[vars[2:]], 99.75) )

       
        
    for var_idx, var in enumerate(vars):

        # print(var_idx, var)

        if len(vars)>1:
            ax = axs[var_idx]
        else:
            ax = axs

        X = df['Time'].values.reshape(size_i, size_j)
        Y = df['Depth'].values.reshape(size_i, size_j)
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
        ##TEST
        #ax.set_xlim(right=20*60)
        ##TEST

        ax.invert_yaxis()

        if var == vars[-1]:
            ax.set_xlabel('time [s]')

        if var == vars[0]:
            secondary_scale(logscale, ax)

    fig.suptitle(title)
    fig.align_ylabels()
    fig.tight_layout(rect=[0,0,1,0.98])
    fig.savefig(r'fig_{:s}.png'.format(title))




def plot_OFT(title, df, items, vars, logscale):
    w = 6
    h = (2.8-1)*len(vars)
    
    w = 8.1
    h = 5.85
    rcParams['figure.figsize'] = [w,h]

    

    fig, axs = plt.subplots(len(vars),1, sharex=True)

    for var_idx, var in enumerate(vars):

        # print(var_idx, var)



        if len(vars)>1:
            ax = axs[var_idx]
        else:
            ax = axs

        for item in items:
            # print('plot item {:d} in {:s} plot'.format(item,var))

            if title =='COFT':
                eleme = ip_mesh['ELEME'].copy()
                eleme = eleme.set_index(('ElName'))

    
                el1 = ip_mesh['CONNE'].loc[item, 'EL1']
                el2 = ip_mesh['CONNE'].loc[item, 'EL2']
                k_dir = ip_mesh['CONNE'].loc[item, 'ISOT']

                mat1 = eleme.loc[el1,'MAT']
                mat2 = eleme.loc[el2,'MAT']
                item_label =  '{:s}>{:s}({:s} to {:s} in {:s} dir.)'.format(el1,el2,mat1, mat2, perm_dict[k_dir])



            elif title =='FOFT':
                el = ip_mesh['ELEME'].loc[item,'ElName']
                mat = ip_mesh['ELEME'].loc[item,'MAT']
                item_label = '{:<4d}{:s}({:s})'.format(item,el,mat)
            df.plot(x='time', y=(item, var), ax=ax, label=item_label, legend = False)
            # df.plot(x='time', y=(item, var), ax=ax, label='tt', legend = False)

        ##TEST
        #ax.set_xlim(left=-1, right=30*60)
        ##TEST

        if var == vars[-1]:
            ax.set_xlabel('time [s]')

        if var == vars[0]:
            secondary_scale(logscale, ax)
        

        
        
        ax.set_ylabel('{:s} [{:s}]'.format(var, units_dict[var]))
        
    if title == 'FOFT':
        ax.set_yscale('log')

    fig.suptitle(title)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right')
    fig.align_ylabels()
    fig.tight_layout(rect=[0,0,0.75,0.98])
    fig.savefig(r'fig_{:s}.png'.format(title))






#Define variables for plotting based on input arguments








if __name__ == '__main__':

    raw_names = [r'FFlow', r'FStatus', r'COFT', r'FOFT']
    fnames = []
    op_files = []
    plot_bool = dict()

    for f in os.listdir(os.getcwd()):
        if f in raw_names:
            f_size = os.path.getsize(f)
            if f_size>0:
                fnames.append(f)
                op_files.append(f)
                plot_bool[f] = True


    args = sys.argv

    ip_file = args[1]




    logscale = False

    if 'log' in args:
        logscale = True
        args.remove('log')



    ffunc = []


    parse_dict = dict()

    for file in fnames:
        if file=='FFlow':
            parse_dict[file]=read_FFlow
        elif file=='FStatus':
            parse_dict[file]=read_FStatus
        elif file=='COFT':
            parse_dict[file]=read_COFT
        elif file=='FOFT':
            parse_dict[file]=read_FOFT

    plot_dict = dict()

    if len(args)>2:
        
        #Turn off all plotting
        for f in plot_bool:
            plot_bool[f] = False


        #Based on arguments turn on the selected files to plot

        for arg in args[2:]:
        


            arg = arg.split('i')




            arg_v = arg[0].strip(",").split(",")

            try:
                arg_i = arg[1].strip(",").split(",")
            except:
                arg_i = None


            
            

            #Check if FStatus is being queried
            if 'FStatus' in fnames  and arg_v[0].lower() == r'FStatus'.lower():
                plot_bool[r'FStatus'] = True
                if len(arg_v)>1:
                    plot_dict[r'FStatus'] = arg_v[1:]
                else:
                    plot_dict[r'FStatus'] = 'all'
                    


            #Check if FFLow is being queried
            elif 'FFlow' in fnames  and arg_v[0].lower() == r'FFlow'.lower():


                plot_bool[r'FFlow'] = True

                if len(arg_v)>1:
                    plot_dict[r'FFlow'] = arg_v[1:]
                else:
                    plot_dict[r'FFlow'] = 'all'



            #Check if COFT is being queried
            elif 'COFT' in fnames  and arg_v[0].lower() == r'COFT'.lower():
                plot_bool[r'COFT'] = True
                plot_dict[r'COFT'] = dict()
                if len(arg_v)>1:
                    plot_dict[r'COFT']['var'] = arg_v[1:]
                else:
                    plot_dict[r'COFT']['var'] = 'all'

                if arg_i is None:
                    plot_dict[r'COFT']['item'] = 'all'
                else:
                    plot_dict[r'COFT']['item'] = arg_i
                

            #Check if FOFT is being queried
            elif 'FOFT' in fnames  and arg_v[0].lower() == r'FOFT'.lower():

                plot_bool[r'FOFT'] = True
                plot_dict[r'FOFT'] = dict()
                if len(arg_v)>1:
                    plot_dict[r'FOFT']['var'] = arg_v[1:]
                else:
                    plot_dict[r'FOFT']['var'] = 'all'

                if arg_i is None:
                    plot_dict[r'FOFT']['item'] = 'all'
                else:
                    plot_dict[r'FOFT']['item'] = arg_i

        


    else:

        for f in plot_bool:

            plot_dict[f] = 'all'

            if f in ['COFT', 'FOFT']:
                plot_dict[f]=dict()
                plot_dict[f]['var'] = 'all'
                plot_dict[f]['item'] = 'all'


    ip_mesh = read_ipMESH(ip_file)

    xl_file = ip_file.split(".")[0]+".xlsx"

    #Query and index data

    for file in fnames:


        plot_f = plot_bool[file]


        if plot_f:
            print('Plotting {:s} data'.format(file))

            if file in ['FFlow', 'FStatus']:
                queried_vars = plot_dict[file]

                df_vars, df = parse_dict[file](file)

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
                
                print('{:s} plot includes: {:s}'.format(file, ' '.join(selected_var)))
                plot_Ffigure(file,df,selected_var)




            elif file in ['COFT', 'FOFT']:
                # print(plot_dict)

                queried_vars = plot_dict[file]['var']
                queried_items = plot_dict[file]['item']
                
                # print(queried_vars)
                
                
                df_vars, df_items, df = parse_dict[file](file)
                

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
                plot_OFT(file, df, selected_items, selected_var)






    spreadsheet = ip_file.split(".")[0]+".xls"
            

    print_Excel = True

    xls_output = input('Do you want to store output files as spreadsheet (Y/N). Default is Yes?\t')

    if xls_output.lower().startswith('n'):
        print_Excel = False
        
    if print_Excel:
        pd_units = pd.Series(units_dict_v2)

        
        if r'FStatus' in fnames:
            #Add column names to FStatus

            fs_var, fs = read_FStatus(r'FStatus')
            fs_row1 = fs.columns.to_list()
            fs_row2 = pd_units[fs_row1].to_list()

            fs_cols = pd.MultiIndex.from_tuples(list(zip( fs_row1, fs_row2)))

            fs.columns = fs_cols

        if r'FFlow' in fnames:

            #Add column names to FFlow

            ff_var, ff = read_FFlow(r'FFlow')
            ff_row1 = ff.columns.to_list()
            ff_row2 = pd_units[ff_row1].to_list()

            ff_cols = pd.MultiIndex.from_tuples(list(zip( ff_row1, ff_row2)))

            ff.columns = ff_cols


        if r'COFT' in fnames:

            #Add column names to COFT
            coft_var, coft_idx, coft = read_COFT(r'COFT')

            eleme2 = ip_mesh['ELEME'].copy()
            eleme2 = eleme2.set_index('ElName')


            query_x = ip_mesh['CONNE']['ISOT']==1
            query_y = ip_mesh['CONNE']['ISOT']==2
            query_z = ip_mesh['CONNE']['ISOT']==3





            coft_row1 = coft.columns.get_level_values(0).to_list()
            coft_row2 = (ip_mesh['CONNE'].loc[coft_row1[1:],'EL1']+ip_mesh['CONNE'].loc[coft_row1[1:],'EL2']).to_list()
            coft_row2 = [''] + coft_row2
            coft_row3 = coft.columns.get_level_values(1).to_list()
            coft_row4 = ['s'] + list(pd_units[coft.columns.get_level_values(1)[1:].to_list()].values)

            coft_cols = pd.MultiIndex.from_tuples(list(zip( coft_row1, 
                                                            coft_row2, 
                                                            coft_row3,
                                                            coft_row4)))

            coft.columns = coft_cols

        if r'FOFT' in fnames:

            #Add column names to FOFT
            foft_var, foft_idx, foft = read_FOFT(r'FOFT')

            foft_row1 = foft.columns.get_level_values(0).to_list()
            foft_row1[0] = 'cell_idx'

            foft_row2 = ip_mesh['ELEME'].loc[foft_row1[1:],'ElName'].to_list()
            foft_row2 = ['cell_name'] + foft_row2

            foft_row2a = ['cell_X [m]'] + ip_mesh['ELEME'].loc[foft_row1[1:], 'X'].to_list()
            foft_row2b = ['cell_Y [m]'] + ip_mesh['ELEME'].loc[foft_row1[1:], 'Y'].to_list()
            foft_row2c = ['cell_Z [m]'] + ip_mesh['ELEME'].loc[foft_row1[1:], 'Z'].to_list()
            foft_row2d = ['cell_mat'] + ip_mesh['ELEME'].loc[foft_row1[1:], 'MAT'].to_list()

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


        with pd.ExcelWriter(xl_file) as writer:
            if r'FStatus' in fnames:
                fs.to_excel(writer, sheet_name=r'FStatus')
            if r'FFlow' in fnames:
                ff.to_excel(writer, sheet_name=r'FFlow')
            if r'COFT' in fnames:
                coft.to_excel(writer, sheet_name=r'COFT')
            if r'FOFT' in fnames:
                foft.to_excel(writer, sheet_name=r'FOFT')


# plt.show()
