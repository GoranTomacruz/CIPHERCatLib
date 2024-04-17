import dearpygui.dearpygui as dpg # type: ignore
import pandas as pd
import numpy as np
import ase
from ase import Atoms
from ase import io
from ase.io import espresso
import numpy as np
import pandas as pd
import os
from os import listdir
from os import mkdir
from os.path import isfile, join, isdir
import os
import ast
from PyAstronomy import pyasl
from ase.io import write

element_list = ['N/A','H','He','Li','Be','B','C','N','O','F','Ne',
           'Na','Mg','Al','Si','P','S','Cl','Ar','K', 'Ca',
           'Sc', 'Ti', 'V','Cr', 'Mn', 'Fe', 'Co', 'Ni',
           'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
           'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
           'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
           'I', 'Xe','Cs', 'Ba','La', 'Ce', 'Pr', 'Nd', 'Pm',
           'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
           'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
           'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
           'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',
           'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
           'Rf', 'Db', 'Sg', 'Bh','Hs', 'Mt', 'Ds', 'Rg', 'Cn',
           'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

occupations_list = ['smearing', 'tetrahedra', 'tetrahedra_lin', 'tetrahedra_opt', 'fixed', 'from_input']
smearing_list = ['gaussian', 'methfessel-paxton', 'marzari-vanderbilt', 'fermi-dirac']
nspin_list = [1, 2, 4]

packing_list = ['fcc', 'bcc', 'hcp']

def get_item_1(sender, app_data, user_data):  
    if app_data != 'N/A':
        dpg.set_value("text_item1", "Searching: {}".format(str(app_data)))  
    else:
        dpg.set_value("text_item1", "Awaiting input on Molecular Search Input...")  
    return app_data

def get_item_2(sender, app_data, user_data):  
    if app_data != 'N/A':  
        dpg.set_value("text_item2", ", {}".format(str(app_data)))  
    else:
        dpg.set_value("text_item2", "")  
    return app_data

def get_item_3(sender, app_data, user_data): 
    if app_data != 'N/A':     
        dpg.set_value("text_item3", ", {}".format(str(app_data)))
    else:
        dpg.set_value("text_item3", "") 
    return app_data   
    
def get_item_4(sender, app_data, user_data):
    if app_data != 'N/A':      
        dpg.set_value("text_item4", ", {}".format(str(app_data)))
    else:
        dpg.set_value("text_item4", "") 
    return app_data   

def get_item_5(sender, app_data, user_data):    
    if app_data != 'N/A':  
        dpg.set_value("text_item5", ", {}".format(str(app_data)))  
    else:
        dpg.set_value("text_item5", "")  
    return app_data

def get_pack(sender, app_data, user_data): 
    if app_data != 'N/A':  
        dpg.set_value("text_item6", "in {}-packing".format(str(app_data)))  
    else:
        dpg.set_value("text_item6", "")   
    return app_data

def _log(sender, app_data, user_data):
        print(f"sender: {sender}, \t app_data: {app_data}, \t user_data: {user_data}")

def clb_selectable(sender, app_data, user_data):
    global compound
    compound = array[user_data][0]
    dpg.delete_item("selected_material")
    dpg.add_text(f"Selecting {array[user_data][0]} in {array[user_data][1]}({array[user_data][2]}). \nPlease make sure input all information \nin the \'Input File Conditions\' window.",tag="selected_material", parent= "input_file_download")

def table_filter(elements,packing):
    df = pd.read_csv('./datasets/Cathub-dataset.csv')                          # Take your df from wherever
    arr = df.to_numpy()                                 ### Convert the DataFrame to a NumPy array.
    i = 0
    while i < len(arr[:]):
        # Test if all elements are present in list
        # Using list comprehension + all()  
        present = all(ele in arr[i][0] for ele in elements)
        if present == 0 or packing not in arr[i][1]:
            arr = np.delete(arr, i, 0)
        else:
            i = i + 1
    return arr

def submit_button():
    elements = []
    if dpg.get_value(elem1) != 'N/A' and len(dpg.get_value(elem1)) != 0 and dpg.get_value(elem1) not in elements:
        elements.append(dpg.get_value(elem1))
    if dpg.get_value(elem2) != 'N/A' and len(dpg.get_value(elem2)) != 0 and dpg.get_value(elem2) not in elements:
        elements.append(dpg.get_value(elem2))
    if dpg.get_value(elem3) != 'N/A' and len(dpg.get_value(elem3)) != 0 and dpg.get_value(elem3) not in elements:
        elements.append(dpg.get_value(elem3))
    if dpg.get_value(elem4) != 'N/A' and len(dpg.get_value(elem4)) != 0 and dpg.get_value(elem4) not in elements:
        elements.append(dpg.get_value(elem4))
    if dpg.get_value(elem5) != 'N/A' and len(dpg.get_value(elem5)) != 0 and dpg.get_value(elem5) not in elements:
        elements.append(dpg.get_value(elem5))
    packing = dpg.get_value(pack)
    
    global array 
    array = table_filter(elements,packing)

    columns = ['Surface','Facet','Miller Indices','Database Source']

    dpg.delete_item("error_statement")
    dpg.delete_item("dataset_table")

    if len(array) == 0:
        dpg.add_text('Search not found, please try again.',tag="error_statement", parent="molecule_results")
    else:
        with dpg.theme() as table_theme:
            with dpg.theme_component(dpg.mvTable):
                dpg.add_theme_color(dpg.mvThemeCol_HeaderActive, (0, 0, 0, 0), category=dpg.mvThemeCat_Core)
                dpg.add_theme_color(dpg.mvThemeCol_Header, (0, 0, 0, 0), category=dpg.mvThemeCat_Core)
        with dpg.table(label='DatasetTable', tag="dataset_table", row_background=True,
                    borders_innerH=True, borders_outerH=True, borders_innerV=True,
                    borders_outerV=True, parent="molecule_results") as material_table:
            for i in range(array.shape[1]):                    # Generates the correct amount of columns
                dpg.add_table_column(label=columns[i])   # Adds the headers
            for i in range(array.shape[0]):                                   # Shows the first n rows
                with dpg.table_row():
                    for j in range(array.shape[1]):
                        dpg.add_selectable(label=f"{array[i,j]}", span_columns=True, callback=clb_selectable, user_data=i)     # Displays the value of
                                                                # each row/column combination
        dpg.bind_item_theme(material_table, table_theme)
    return array

def submit_button_2():
    print(compound)
    calc_style = dpg.get_value("calc_style")
    ecutwfc = dpg.get_value("ecutwfc")
    occupations = dpg.get_value("occupations")
    smearing = dpg.get_value("smearing")
    nspin = dpg.get_value("nspin")
    kpt_x = dpg.get_value("kpt_x")
    kpt_y = dpg.get_value("kpt_y")
    kpt_z = dpg.get_value("kpt_z")
    if '6' in compound:
        compound_2 = compound.replace('6','2')
        sym_1 = compound_2[0:find_nth(compound_2,'2',1)]
        sym_2 = compound_2[find_nth(compound_2,'2',1)+1:find_nth(compound_2,'2',2)]
        compound_3 = sym_2 + '2' + sym_1 + '2'
    else:
        compound_2 = compound
        compound_3 = compound
    make_input_file(compound_2,compound_3,calc_style,ecutwfc,occupations,smearing,nspin,kpt_x,kpt_y,kpt_z)

#FUNCTIONS FOR REFINING INPUT FILE
def replacer(symbol,filename):
    f = open(filename,'r')
    filedata = f.read()
    f.close()
    newdata = filedata.replace("\'0.01D0\'","0.01D0")
    newdata = newdata.replace("\'0.0D0\'","0.0D0")
    newdata = newdata.replace("'.true.'",".true.")
    f = open(filename,'w')
    f.write(newdata)
    f.close()

def find_nth(string,substring,n):
    if n == 1:
        return string.find(substring)
    else:
        return string.find(substring, find_nth(string, substring, n - 1) + 1)
    
# INPUT FILE MAKING
def make_input_file(compound,alt_compound,calc_style,ecutwfc,occupations,smearing,nspin,kpt_x,kpt_y,kpt_z):
    filecounter = 1

    an = pyasl.AtomicNo()
    if isdir('./input_files') == False:
        mkdir('./input_files')
    if isdir('./pictures') == False:
        mkdir('./pictures')

    onlyfiles = [f for f in listdir('structures/') if isfile(join('structures/', f))]
    finalfiles = []

    text_option_1 = "#" + compound + '.txt'
    text_option_2 = "_" + compound + '.txt'
    text_option_3 = "_" + alt_compound + '.txt'

    # filter to get only the output files
    for file in onlyfiles:
        if text_option_1 == file:
            finalfiles.append(file)
        elif text_option_2 == file or text_option_3 == file:
            finalfiles.append(file)

    errorfiles = []
    exceptions = []
    coords_arr = []
    axes_arr = []
    finalfiles = sorted(finalfiles)
    print(finalfiles)

    for file in finalfiles:
        sym = file[1:file.find('.t')]
        print(sym)
        file = 'structures/' + file
        fileScan = open(file, "r")
        
        x = 0
        y = 0
        z = 0
        coord = []
        axes = []
        symbols = []
        crystal = 0

        for line in fileScan:
            if 'Coordinates:' in line or (x > 0 and x <= 12):    
                if x > 0:
                    coord.append([float(i) for i in line.split()])
                x = x + 1

            if 'Crystal Axes:' in line or (y > 0 and y <= 3):    
                if y > 0:
                    axes.append([float(i) for i in line.split()])
                y = y + 1
            if 'Order:' in line or (z > 0 and z <= 1):    
                if z > 0:
                    symbols = ast.literal_eval(line)
                    symbols = [n.strip() for n in symbols]
                z = z + 1

        if symbols == []:
            symbols = [sym]*12

        print(symbols)

        coords_arr = np.array(coord)
        axes_arr = np.array(axes)

        #if any coordinate is greater than 1, crystal = 1
        if (coords_arr < 1.05).all():
            crystal = 1

        if crystal == 1:
            #print('conversion')
            for i in range(len(coords_arr)):
                coords_arr[i] = np.matmul(coords_arr[i],axes_arr)

        coords_arr = list(zip(coords_arr.T[0],coords_arr.T[1],coords_arr.T[2]))
        axes_arr = list(zip(axes_arr.T[0],axes_arr.T[1],axes_arr.T[2]))

        filecounter = filecounter + 1
        ecutwfc_ry = ecutwfc*0.0734986176495
        ecutrho_ry = ecutwfc_ry*10

        filename = './input_files/' + sym + "-" + calc_style + ".in"
        picname = './pictures/' + sym + ".png"
        outfilename = './input_files/' + sym + ".out"
        symbol_count = len(symbols)
        numbers = np.zeros([12], dtype=int)
        kpoints = (kpt_x, kpt_y, kpt_z)

        i = 1
        while i <= len(symbols):
            numbers[i-1] = an.getAtomicNo(symbols[i-1])
            i = i + 1

        pseudo = {
        "Sc": "Sc.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Ti": "Ti.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "V": "V.pbe-spnl-rrkjus_psl.1.0.0.UPF",
        "Cr": "Cr.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Mn": "Mn.pbe-spn-rrkjus_psl.0.3.1.UPF",
        "Fe": "Fe.pbe-n-rrkjus_psl.1.0.0.UPF",
        "Co": "Co.pbe-spn-rrkjus_psl.0.3.1.UPF",
        "Ni": "Ni.pbe-n-rrkjus_psl.1.0.0.UPF",
        "Cu": "Cu.pbe-dn-rrkjus_psl.1.0.0.UPF",
        "Zn": "Zn.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Y": "Y.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Zr": "Zr.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Nb": "Nb.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Mo": "Mo.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Tc": "Tc.pbe-spn-rrkjus_psl.0.3.0.UPF",
        "Ru": "Ru.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Rh": "Rh.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Pd": "Pd.pbe-n-rrkjus_psl.1.0.0.UPF",
        "Ag": "Ag.pbe-n-rrkjus_psl.1.0.0.UPF",
        "Cd": "Cd.pbe-n-rrkjus_psl.1.0.0.UPF",
        "La": "La.pbe-spfn-rrkjus_psl.1.0.0.UPF",
        "Hf": "Hf.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Ta": "Ta.pbe-spfn-rrkjus_psl.1.0.0.UPF",
        "W": "W.pbe-spn-rrkjus_psl.0.2.3.UPF",
        "Re": "Re.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Os": "Os.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Ir": "Ir.pbe-n-rrkjus_psl.0.2.3.UPF",
        "Pt": "Pt.pbe-n-rrkjus_psl.1.0.0.UPF",
        "Au": "Au.pbe-n-rrkjus_psl.1.0.0.UPF",
        "Hg": "Hg.pbe-n-rrkjus_psl.1.0.0.UPF",
        "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
        "O": "O.pbe-n-rrkjus_psl.1.0.0.UPF",
        "C": "C.pbe-n-rrkjus_psl.1.0.0.UPF"
        }

        structure = Atoms(numbers, positions=coords_arr, cell=axes_arr)

        if calc_style == 'scf' or calc_style == 'nscf':
            inputdata = {
            'calculation':calc_style,
            'tstress':'.true.',
            'tprnfor':'.true.',
            'tefield':'.true.',
            'dipfield':'.true.',
            'pseudo_dir':'../pseudo',
            'outdir':'./' + sym,
            'prefix': sym + '_beef',
            'ecutwfc':ecutwfc_ry,
            'ecutrho':ecutrho_ry,
            'nosym':'.true.',
            'occupations': occupations,
            'smearing': smearing,
            'nspin': int(nspin),
            'input_dft':'BEEF-vdW',
            'edir':3,
            'emaxpos':0.95,
            'eopreg':0.05,
            'eamp': 0.01,
            'conv_thr':0.00001,
            'mixing_beta':0.1,
            'electron_maxstep':500,
            'mixing_mode':'plain',
            'mixing_ndim':10,
            'degauss':'0.01D0',
            'celldm(1)':1.889726
            }
        elif calc_style == 'relax':
            inputdata = {
            'calculation':calc_style,
            'ion_dynamics':'bfgs',
            'tstress':'.true.',
            'tprnfor':'.true.',
            'tefield':'.true.',
            'dipfield':'.true.',
            'pseudo_dir':'../pseudo',
            'outdir':'./' + sym,
            'prefix': sym + '_beef',
            'ecutwfc':ecutwfc_ry,
            'ecutrho':ecutrho_ry,
            'nosym':'.true.',
            'occupations': occupations,
            'smearing': smearing,
            'nspin': int(nspin),
            'input_dft':'BEEF-vdW',
            'edir':3,
            'emaxpos':0.95,
            'eopreg':0.05,
            'eamp': 0.01,
            'conv_thr':0.00001,
            'mixing_beta':0.1,
            'electron_maxstep':500,
            'mixing_mode':'plain',
            'mixing_ndim':10,
            'degauss':'0.01D0',
            'celldm(1)':1.889726
            }

        f = open(filename, "w")

        ase.io.espresso.write_espresso_in(f, structure, input_data=inputdata, pseudopotentials=pseudo, kspacing=None, kpts=kpoints, koffset=(0,0,0))
        f.close()

        write(picname, structure)

        replacer(sym,filename)
        dpg.delete_item("selected_material")
        dpg.delete_item("molecule_structure")
        dpg.delete_item("input_file_text")
        dpg.add_text(f"Unit Cell Figure of {compound}:",tag="selected_material", parent= "input_file_download")
        add_and_load_image(picname, parent='input_file_download')
        dpg.add_text(f"The input File of {compound} is in your folder! \nSimply go to {filename} to access it.\nBest of luck in your runs!",tag="input_file_text", parent= "input_file_download", pos=(10,250))
       
#ADDITION OF PICTURE
def add_and_load_image(image_path, parent):
    width, height, channels, data = dpg.load_image(image_path)

    with dpg.texture_registry() as reg_id:
        texture_id = dpg.add_static_texture(width, height, data, parent=reg_id)

    if parent is None:
        return dpg.add_image(texture_id, tag='molecule_structure', pos=(150,50))
    else:
        return dpg.add_image(texture_id, parent=parent, tag='molecule_structure', pos=(150,50))


#START OF DEARPYGUI MAKING
dpg.create_context()



with dpg.window(label="Molecule Search User Input", pos=(0,0), width=500, height=400, no_collapse=True, tag='molecular_search_user_input'):
    dpg.add_text("Enter elements to search for:")
    elem1 = dpg.add_combo(label="Element 1", items = element_list, callback=get_item_1)
    elem2 = dpg.add_combo(label="Element 2", items = element_list, callback=get_item_2)
    elem3 = dpg.add_combo(label="Element 3", items = element_list, callback=get_item_3)
    elem4 = dpg.add_combo(label="Element 4", items = element_list, callback=get_item_4)
    elem5 = dpg.add_combo(label="Element 5", items = element_list, callback=get_item_5)
    pack = dpg.add_combo(label="Select Atomic Packing", items = packing_list, callback=get_pack) 
    dpg.add_button(label="Submit!", callback=submit_button)
    pass

with dpg.window(label="Molecule Results", tag="molecule_results", pos=(500,0), width=500, height=400, no_collapse=True):
    with dpg.group(horizontal=True):
        dpg.add_text("Awaiting input on Molecular Search Input... ", tag="text_item1")
        dpg.add_text("", tag="text_item2")
        dpg.add_text("", tag="text_item3")
        dpg.add_text("", tag="text_item4")
        dpg.add_text("", tag="text_item5")
        dpg.add_text("", tag="text_item6")
    pass

with dpg.window(label="Input File Conditions", pos=(0,400), width=500, height=400, no_collapse=True, tag='input_file_conditions'):
    dpg.add_text("QE Input")
    with dpg.group(horizontal=True):
        dpg.add_text("Calculation: ")
        dpg.add_combo(items = ['relax', 'scf', 'nscf'], tag='calc_style')
    with dpg.group(horizontal=True):
        dpg.add_text("Cutoff energy for wavefunction (eV): ")
        dpg.add_input_float(width = 250, callback=_log, step=0, tag='ecutwfc')
    pass
    with dpg.group(horizontal=True):
        dpg.add_text("Occupations: ")
        dpg.add_combo(items = occupations_list, tag='occupations')
    pass
    with dpg.group(horizontal=True):
        dpg.add_text("Smearing: ")
        dpg.add_combo(items = smearing_list, tag='smearing')
    with dpg.group(horizontal=True):
        dpg.add_text("n-spin: ")
        dpg.add_combo(items = nspin_list, tag='nspin')
    with dpg.group(horizontal=True):
        dpg.add_text("k-points: ")
        dpg.add_input_int(callback=_log, width = 30, step=0, tag='kpt_x')
        dpg.add_text(" x ")
        dpg.add_input_int(callback=_log, width = 30, step=0, tag='kpt_y')
        dpg.add_text(" x ")
        dpg.add_input_int(callback=_log, width = 30, step=0, tag='kpt_z')
    dpg.add_button(label="Submit!", callback=submit_button_2)
    pass

with dpg.window(label="Input File Download", pos=(500,400), width=500, height=400, no_collapse=True, tag = "input_file_download"):
    pass

dpg.create_viewport(title='CatLib Demo Alpha Version', width=1000, height=800)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()