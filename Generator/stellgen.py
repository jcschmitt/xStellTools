import tkinter as tk
from tkinter import ttk, filedialog, messagebox
#import tkinter.filedialog as fd
import os
import shutil
#import sys
import f90nml
import numpy as np
from numpy import iinfo


class stellgen:
# 'Notebook'-style gui
# Front frame is main control frame
#       Export, Import, etc.
#       Output name is a function of backframe settings
# Front frame has loop control
# Back frames for VMEC, REGCOIL, OPTIMUM, etc
#    Pre-made settings should be avail on each frame
#    

    def __init__(self, my_root):
        # Set up the main window, its frames, and their components
        self.my_root = my_root

        #########################
        # set up the tabs(menus)
        #########################

        # For part 1 (Graphical User Interface controls and support functions)
        #Background colorscheme
        self.bg_color_1 = "#605d20"
        self.bg_color_2 = "#cc8888"
        self.bg_color_3 = "#b0c4de"
        
        self.tab_control = ttk.Notebook(self.my_root)
        self.vmec_tab = ttk.Frame(self.tab_control)
        self.optimum_tab = ttk.Frame(self.tab_control)
        self.targets_tab = ttk.Frame(self.tab_control)
        self.variables_tab = ttk.Frame(self.tab_control)
        self.regcoil_tab = ttk.Frame(self.tab_control)
        self.bootstrap_tab = ttk.Frame(self.tab_control)
        self.scanables_tab = ttk.Frame(self.tab_control)
        self.filenames_tab = ttk.Frame(self.tab_control)

        
        #self.make_main_tab(self.main_tab)
        self.make_filenames_tab(self.filenames_tab)
        self.make_optimum_tab(self.optimum_tab)
        self.make_targets_tab(self.targets_tab)
        self.make_variables_tab(self.variables_tab)
        self.make_scanables_tab(self.scanables_tab)
        self.make_vmec_tab(self.vmec_tab)
        self.make_regcoil_tab(self.regcoil_tab)
        self.make_bootstrap_tab(self.bootstrap_tab)

        self.tab_control.add(self.filenames_tab, text='Files')
        self.tab_control.add(self.scanables_tab, text='Scan Control')
        self.tab_control.add(self.optimum_tab, text='Optimum')
        self.tab_control.add(self.targets_tab, text='Targets')
        self.tab_control.add(self.variables_tab, text='Variables')
        self.tab_control.add(self.vmec_tab,
                            text='VMEC')
        self.tab_control.add(self.regcoil_tab, text='Regcoil')
        self.tab_control.add(self.bootstrap_tab, text='Profiles &\nBootstrap')
        self.tab_control.grid()

        my_root.title("STELLOPT Generator 0.1a")

        the_top = my_root.winfo_toplevel()
        the_menu_bar = tk.Menu(the_top)
        the_top['menu'] = the_menu_bar


    def add_label_and_entry_to_frame(self, position,
                                     in_key, in_dict, target_frame):
        entry_width = 12
        if (type(in_dict[in_key]) is tk.StringVar):
            entry_width = max(entry_width, len(in_dict[in_key].get()) + 4)

        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.E,
                      text=(in_key+":"))
        this_lbl.grid(row=position, rowspan=1, column=1, columnspan=1,
                      sticky=tk.E)

        this_entry = tk.Entry(target_frame,
                              width=entry_width,
                              textvariable=in_dict[in_key])
        this_entry.grid(row=position, rowspan=1, column=2, columnspan=1,
                        sticky=tk.W)

    def add_label_and_scalar_param_to_frame(self, position,
                                     in_key, in_dict, target_frame):
        entry_width = 12

        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.E,
                      text=(in_key+":"))
        this_lbl.grid(row=position, rowspan=1, column=1, columnspan=1,
                      sticky=tk.E)

        this_entry = tk.Entry(target_frame,
                              width=entry_width,
                              textvariable=in_dict[in_key])
        this_entry.grid(row=position, rowspan=1, column=2, columnspan=1,
                        sticky=tk.W)
   
    def add_text_entry_to_frame(self, position,
                                in_key, in_dict, target_frame,
                                entry_width=30, entry_height=14):

        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.CENTER,
                      text=(in_key+":"))
        this_lbl.grid(row=position, rowspan=1, column=1, columnspan=1,
                      sticky=tk.E)

        this_entry = tk.Text(target_frame,
                             height=entry_height,
                             width=entry_width)

        this_entry.grid(row=position, rowspan=1, column=2, columnspan=1,
                        sticky=tk.W)

        this_entry.insert('1.0', in_dict[in_key].get())

        return this_entry

    def new_text_entry_data(self, position,
                            in_key, in_dict, target_frame,
                            entry_width=30, entry_height=14):

        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.CENTER,
                      text=(in_key+":"))
        this_lbl.grid(row=position, rowspan=1, column=1, columnspan=1,
                      sticky=tk.E)

        this_entry = tk.Text(target_frame,
                             height=entry_height,
                             width=entry_width)

        this_entry.grid(row=position, rowspan=1, column=2, columnspan=1,
                        sticky=tk.W)

        this_entry.insert('1.0', in_dict[in_key].get())

        return this_entry

    def add_label_and_scan_line_to_frame(self, position,
                                     in_key, in_dict, target_frame, with_weights=False,
                                     entry_width=12, col_span=1):
    

#         this_lbl = tk.Label(target_frame,
#                       bg=self.bg_color_1,
#                       anchor=tk.E,
#                       text=(in_key+":"))
#         this_lbl.grid(row=position, rowspan=1, column=1, columnspan=1,
#                       sticky=tk.E)

        this_rb = tk.Checkbutton(target_frame,
                                       bg=self.bg_color_1,
                                       text=in_key,
                                       variable=in_dict[in_key]['Enabled'])
        this_rb.grid(row=position, rowspan=1, column=1, columnspan=1, sticky=tk.W)

        this_entry = tk.Entry(target_frame,
                              width=entry_width,
                              textvariable=in_dict[in_key]['Values'])
        this_entry.grid(row=position, rowspan=1, column=2, columnspan=col_span,
                        sticky=tk.W)

        if (with_weights is True):
            this_entry = tk.Entry(target_frame,
                                  width=entry_width,
                                  textvariable=in_dict[in_key]['Sigmas'])
            this_entry.grid(row=position, rowspan=1, column=2+col_span, columnspan=col_span,
                            sticky=tk.W)

    def add_selectables_line_to_frame(self, position1, position2,
                                     in_key, in_dict, target_frame):

        this_rb = tk.Checkbutton(target_frame,
                                       bg=self.bg_color_1,
                                       text=in_key,
                                       variable=in_dict[in_key]['Enabled'])
        this_rb.grid(row=position1, rowspan=1, column=position2, columnspan=1, sticky=tk.W)




    def add_optimum_legend_to_frame(self, position, target_frame):
        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.E,
                      text='Target')
        this_lbl.grid(row=position, rowspan=1, column=2, columnspan=1,
                      sticky=tk.E)
        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.E,
                      text='Sigma')
        this_lbl.grid(row=position, rowspan=1, column=3, columnspan=1,
                      sticky=tk.E)
        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.E,
                      text='Count')
        this_lbl.grid(row=position, rowspan=1, column=4, columnspan=1,
                      sticky=tk.E)
        
    

    def add_optimum_target_entry_to_frame(self, position,
                                     in_key, in_dict, target_frame):
        entry_width = 12

        this_lbl = tk.Label(target_frame,
                      bg=self.bg_color_1,
                      anchor=tk.E,
                      text=(in_key+":"))
        this_lbl.grid(row=position, rowspan=1, column=1, columnspan=1,
                      sticky=tk.E)

        this_target = tk.Entry(target_frame,
                              width=entry_width,
                              textvariable=in_dict[in_key]['TARGET'])
        this_target.grid(row=position, rowspan=1, column=2, columnspan=1,
                        sticky=tk.W)

        this_sigma = tk.Entry(target_frame,
                              width=entry_width,
                              textvariable=in_dict[in_key]['SIGMA'])
        this_sigma.grid(row=position, rowspan=1, column=3, columnspan=1,
                        sticky=tk.W)

        this_count = tk.Entry(target_frame,
                              width=entry_width,
                              textvariable=in_dict[in_key]['COUNT'])
        this_count.grid(row=position, rowspan=1, column=4, columnspan=1,
                        sticky=tk.W)
        
   
    def make_vmec_tab(self, this_tab):
        # VMEC Execution parameters.
        # Store everything in a dictionary for later
        self.VMEC_RUN_PARAMS = {}
        self.VMEC_RUN_PARAMS['DELT'] = tk.DoubleVar()
        self.VMEC_RUN_PARAMS['DELT'].set(0.9)

        self.VMEC_RUN_PARAMS['TCON0'] = tk.DoubleVar()
        self.VMEC_RUN_PARAMS['TCON0'].set(2.0)

        self.VMEC_RUN_PARAMS['NS_ARRAY'] = tk.StringVar()
        self.VMEC_RUN_PARAMS['NS_ARRAY'].set('51')

        self.VMEC_RUN_PARAMS['FTOL_ARRAY'] = tk.StringVar()
        self.VMEC_RUN_PARAMS['FTOL_ARRAY'].set('1.0E-19')

        self.VMEC_RUN_PARAMS['NTOR'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['NTOR'].set(8)

        self.VMEC_RUN_PARAMS['MPOL'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['MPOL'].set(8)

        self.VMEC_RUN_PARAMS['NITER_ARRAY'] = tk.StringVar()
        self.VMEC_RUN_PARAMS['NITER_ARRAY'].set('10000')

        self.VMEC_RUN_PARAMS['NSTEP'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['NSTEP'].set(200)

        self.VMEC_RUN_PARAMS['NVACSKIP'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['NVACSKIP'].set(10)

        self.VMEC_RUN_PARAMS['LFORBAL'] = tk.BooleanVar()
        self.VMEC_RUN_PARAMS['LFORBAL'].set(False)

        self.VMEC_RUN_PARAMS['LASYM'] = tk.BooleanVar()
        self.VMEC_RUN_PARAMS['LASYM'].set(False)

        self.VMEC_RUN_PARAMS['NZETA'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['NZETA'].set(36)

        self.VMEC_RUN_PARAMS['NTHETA'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['NTHETA'].set([])

        self.VMEC_RUN_PARAMS['NFP'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['NFP'].set(5)

        self.VMEC_RUN_PARAMS['PHIEDGE'] = tk.DoubleVar()
        self.VMEC_RUN_PARAMS['PHIEDGE'].set(0.1)

        self.VMEC_RUN_PARAMS['NCURR'] = tk.IntVar()
        self.VMEC_RUN_PARAMS['NCURR'].set(1)

        self.VMEC_RUN_PARAMS['LFREEB'] = tk.BooleanVar()
        self.VMEC_RUN_PARAMS['LFREEB'].set(False)

        self.VMEC_RUN_PARAMS['MGRID_FILE'] = tk.StringVar()
        self.VMEC_RUN_PARAMS['MGRID_FILE'].set('')

        # needs to be fixed to handle arrays
        self.VMEC_RUN_PARAMS['EXTCUR'] = tk.StringVar()
        self.VMEC_RUN_PARAMS['EXTCUR'].set('1')

        # These parameters need to be written out later.
        self.VMEC_RUN_PARAMS['InitPosition'] = tk.StringVar()
        init_position_str = ('! Initial Position.\n' +
                              '  RAXIS_CC = 5.5607E+00\n' +
                              '  ZAXIS_CS = 0.0000E+00\n' +
                              '  RBC(0,0) = 5.5210E+00\n' +
                              '  RBC(1,0) = 2.7849e-01\n' +
                              '!  RBC(5,0) = 2.7849e-01\n' +
                              '  RBC(0,1) = 4.8900E-01\n' +
                              '  ZBS(0,0) = 0.0000E+00\n' +
                              '  ZBS(1,0) = -2.3504e-01\n' +
                              '!  ZBS(5,0) = -2.3504e-01\n' +
                              '  ZBS(0,1) = 6.2496E-01\n')
        self.VMEC_RUN_PARAMS['InitPosition'].set(init_position_str)

        # do stuff and things
        # Make two frames. 
        vmec_frame1 = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="VMEC Execution Parameters",
                                   font=('Helvetica', '14'))
        vmec_frame1.grid(row=1,
                        rowspan=2,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        vmec_frame2 = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="VMEC Boundary and Axis Parameters",
                                   font=('Helvetica', '14'))
        vmec_frame2.grid(row=1,
                        rowspan=2,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        self.vmec_frame3 = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="VMEC Templates",
                                   font=('Helvetica', '14'))
        self.vmec_frame3.grid(row=3,
                        rowspan=2,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        # Loop over all of the VMEC items and add them to the VMEC frame
        # vmec_keys = self.VMEC_RUN_PARAMS.keys()
        # Not looping because the order changes at time. Listing the explicetly
        counter = 0
        for this_key in ('NTOR', 'MPOL', 'NZETA', 'NTHETA', 'NS_ARRAY', 'FTOL_ARRAY',
                         'NITER_ARRAY', 'NCURR', 'NSTEP', 'NVACSKIP', 'DELT', 'PHIEDGE',
                         'TCON0', 'LFORBAL', 'LASYM', 'NFP', 'LFREEB', 'MGRID_FILE',
                         'EXTCUR'):
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.VMEC_RUN_PARAMS,
                                              vmec_frame1)

        row_counter = 0 # frame 2 will be used
        self.VMEC_InitPosition = self.add_text_entry_to_frame(row_counter,
                                    'InitPosition',
                                    self.VMEC_RUN_PARAMS,
                                    vmec_frame2, entry_width=120)


        row_counter = 1 # frame 3 will be used
        col_counter = 1
        vmec_aten_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_aten, # self.doit,
                                text='ATEN')
        vmec_aten_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_mljs2_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_mljs2, # self.doit,
                                text='MLJS2')
        vmec_mljs2_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_mljs3_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_mljs3, # self.doit,
                                text='MLJS3')
        vmec_mljs3_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_mljs4_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_mljs4, # self.doit,
                                text='MLJS4')
        vmec_mljs4_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_ware5p2_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_ware5p2, # self.doit,
                                text='WARE5P2')
        vmec_ware5p2_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_wb24_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_wb24, # self.doit,
                                text='WB24')
        vmec_wb24_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        row_counter += 1 # frame 3 will be used
        col_counter = 1
        vmec_aten_b1_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_aten_beta_1p7, # self.doit,
                                text='ATEN 1.7%')
        vmec_aten_b1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_mljs2_b1_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_mljs2_beta_1p8, # self.doit,
                                text='MLJS2 1.8%')
        vmec_mljs2_b1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_mljs3_b1_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_mljs3_beta_1p8, # self.doit,
                                text='MLJS3 1.8%')
        vmec_mljs3_b1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_mljs4_b1_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_mljs4_beta_2p1, # self.doit,
                                text='MLJS4 2.1%')
        vmec_mljs4_b1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_ware5p2_b1_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_ware5p2_beta_1p9, # self.doit,
                                text='WARE5P2 1.9%')
        vmec_ware5p2_b1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        vmec_wb24_b1_button = tk.Button(self.vmec_frame3,
                                command=self.load_vmec_init_position_wb24_beta_1p6, # self.doit,
                                text='WB24 1.6%')
        vmec_wb24_b1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

    def load_vmec_p1(self, filename):
        vmec_nml = f90nml.read(filename)
        # numeric values
        for variable in ('NTOR', 'MPOL', 'NZETA', 'NTHETA', 
                         'NCURR', 'NSTEP', 'NVACSKIP', 'DELT',
                         'TCON0', 'NFP', 'MGRID_FILE', 'PHIEDGE'):
            print('<---Loading: ' + variable)
            try:
                in_val = vmec_nml['indata'][variable]
            except:
                print('<----: ' + variable + ' not found')
                in_val = []
            self.VMEC_RUN_PARAMS[variable].set(in_val)

        for variable in ('NS_ARRAY', 'FTOL_ARRAY', 'NITER_ARRAY', 'EXTCUR'):
            print('<---Loading: ' + variable)
            try:
                in_val = vmec_nml['indata'][variable]
                new_str = str(in_val)
                new_str = new_str.replace('[', '')
                new_str = new_str.replace(']', '')
            except:
                print('<----: ' + variable + ' not found')
                in_val = []
                new_str = ''
            self.VMEC_RUN_PARAMS[variable].set(new_str)
                            
        for variable in ('LFORBAL', 'LASYM'):
            print('<---Loading: ' + variable)
            try:
                in_val = vmec_nml['indata'][variable]
            except:
                print('<----: ' + variable + ' not found')
                in_val = False
            self.VMEC_RUN_PARAMS[variable].set(in_val)

            
        #for variable in ('NS_ARRAY', 'FTOL_ARRAY'):
        #    in_val = vmec_nml['indata'][variable]
        #    if ((type(in_val) == int) or (type(in_val) == float)):
        #        self.VMEC_RUN_PARAMS[variable].set(in_val)
        #    else:
        #        self.VMEC_RUN_PARAMS[variable].set(in_val[0])

        

    def load_vmec_p2(self, filename):
        file_handle = open(filename);
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.VMEC_InitPosition.delete('1.0', 'end')
        self.VMEC_InitPosition.insert('1.0', new_text.get())


    def load_vmec_init_position_aten(self):
        self.load_vmec_p1('templates/vmec_aten_a3b25_p1')
        self.load_vmec_p2('templates/vmec_aten_a3b25_p2')
        

    def load_vmec_init_position_aten_beta_1p7(self):
        self.load_vmec_p1('templates/vmec_aten_beta_1p7_p1')
        self.load_vmec_p2('templates/vmec_aten_beta_1p7_p2')
        

    def load_vmec_init_position_mljs2(self):
        self.load_vmec_p1('templates/vmec_mljs2_a3b25_p1')
        self.load_vmec_p2('templates/vmec_mljs2_a3b25_p2')


    def load_vmec_init_position_mljs2_beta_1p8(self):
        self.load_vmec_p1('templates/vmec_mljs2_beta_1p8_p1')
        self.load_vmec_p2('templates/vmec_mljs2_beta_1p8_p2')


    def load_vmec_init_position_mljs3(self):
        self.load_vmec_p1('templates/vmec_mljs3_a3b25_p1')
        self.load_vmec_p2('templates/vmec_mljs3_a3b25_p2')


    def load_vmec_init_position_mljs3_beta_1p8(self):
        self.load_vmec_p1('templates/vmec_mljs3_beta_1p8_p1')
        self.load_vmec_p2('templates/vmec_mljs3_beta_1p8_p2')


    def load_vmec_init_position_mljs4(self):
        self.load_vmec_p1('templates/vmec_mljs4_a3b25_p1')
        self.load_vmec_p2('templates/vmec_mljs4_a3b25_p2')


    def load_vmec_init_position_mljs4_beta_2p1(self):
        self.load_vmec_p1('templates/vmec_mljs4_beta_2p1_p1')
        self.load_vmec_p2('templates/vmec_mljs4_beta_2p1_p2')


    def load_vmec_init_position_ware5p2(self):
        self.load_vmec_p1('templates/vmec_ware5p2_a3b25_p1')
        self.load_vmec_p2('templates/vmec_ware5p2_a3b25_p2')


    def load_vmec_init_position_ware5p2_beta_1p9(self):
        self.load_vmec_p1('templates/vmec_ware5p2_beta_1p9_p1')
        self.load_vmec_p2('templates/vmec_ware5p2_beta_1p9_p2')


    def load_vmec_init_position_wb24(self):
        self.load_vmec_p1('templates/vmec_wb24_a3b25_p1')
        self.load_vmec_p2('templates/vmec_wb24_a3b25_p2')

   
    def load_vmec_init_position_wb24_beta_1p6(self):
        self.load_vmec_p1('templates/vmec_wb24_beta_1p6_p1')
        self.load_vmec_p2('templates/vmec_wb24_beta_1p6_p2')

   
    def make_optimum_tab(self, this_tab):
        # Store everything in a dictionary for later
        self.OPTIMUM_PARAMS = {}
        self.OPTIMUM_PARAMS['NFUNC_MAX'] = tk.IntVar()
        self.OPTIMUM_PARAMS['NFUNC_MAX'].set(5000)

        self.OPTIMUM_PARAMS['EQUIL_TYPE'] = tk.StringVar()
        self.OPTIMUM_PARAMS['EQUIL_TYPE'].set('VMEC2000_ONEEQ')

        self.OPTIMUM_PARAMS['BOOTCALC_TYPE'] = tk.StringVar()
        self.OPTIMUM_PARAMS['BOOTCALC_TYPE'].set('BOOTSJ')

        self.OPTIMUM_PARAMS['OPT_TYPE'] = tk.StringVar()
        self.OPTIMUM_PARAMS['OPT_TYPE'].set('LMDIF_bounded')

        self.OPTIMUM_PARAMS['FTOL'] = tk.DoubleVar()
        self.OPTIMUM_PARAMS['FTOL'].set(1.00e-6)

        self.OPTIMUM_PARAMS['XTOL'] = tk.DoubleVar()
        self.OPTIMUM_PARAMS['XTOL'].set(1.00e-6)

        self.OPTIMUM_PARAMS['GTOL'] = tk.DoubleVar()
        self.OPTIMUM_PARAMS['GTOL'].set(1.00e-30)

        self.OPTIMUM_PARAMS['FACTOR'] = tk.DoubleVar()
        self.OPTIMUM_PARAMS['FACTOR'].set(10.0)

        self.OPTIMUM_PARAMS['EPSFCN'] = tk.DoubleVar()
        self.OPTIMUM_PARAMS['EPSFCN'].set(1.00e-6)

        self.OPTIMUM_PARAMS['MODE'] = tk.IntVar()
        self.OPTIMUM_PARAMS['MODE'].set(1)

        self.OPTIMUM_PARAMS['LKEEP_MINS'] = tk.BooleanVar()
        self.OPTIMUM_PARAMS['LKEEP_MINS'].set(True)

       # These parameters need to be written out later.
        self.OPTIMUM_PARAMS['Extra_Lines'] = tk.StringVar()
        extra_lines_str = (' ! Extra Lines Go Here\n' +
                           ' !  M_LBFGSB = 5\n' +
                           ' !  PRINT_LBFGSB = 1\n' +
                           " ! For Fourier-series winding surface\n" +
                           " ! This is the initial winding surface\n" +
                           "  REGCOIL_NESCIN_FILENAME= 'nescin.out' \n" + 
                           "  REGCOIL_NUM_FIELD_PERIODS = 5\n\n")
        self.OPTIMUM_PARAMS['Extra_Lines'].set(extra_lines_str)

       # These parameters need to be written out later.
        self.OPTIMUM_PARAMS['Profile_Functions'] = tk.StringVar()
        extra_lines_str = ("!---------------------------\n" + 
                           "!       Profile Functions   \n" + 
                           "!---------------------------\n\n" + 
                           "  ! Note that ne_opt is normalized to 1e18 meters^{-3}.\n" + 
                           "  ! n = (0.7e20/m^3) * (1 - s^5)\n" + 
                           "  NE_TYPE = 'power_series'\n" + 
                           "  NE_OPT = 70.0 0.0 0.0 0.0 0.0 -70.0\n\n" + 
                           "! TE_OPT and TI_OPT are in units of 1 eV.\n" + 
                           " ! T  = 2 keV * (1 - s)\n" + 
                           "  TE_TYPE = 'power_series'\n" + 
                           "  TE_OPT = 2e3 -2e3\n" + 
                           "  TI_TYPE = 'power_series'\n" + 
                           "  TI_OPT = 2e3 -2e3\n\n" + 
                           "bootj_type='power_series'\n" + 
                           "! The number of nonzero entries in bootj_aux_f sets the degree of the polynomial fit!\n" + 
                           "bootj_aux_f = 16*1.0e-10\n\n" + 
                           "sfincs_s = 0.00851345, 0.0337639, 0.0748914, 0.130496, 0.198683, 0.277131, 0.363169, 0.453866, 0.546134, 0.636831, 0.722869, 0.801317, 0.869504, 0.925109, 0.966236, 0.991487 ! 16 points\n\n" + 
                           " sfincs_min_procs = 32\n" + 
                           " vboot_tolerance = 1.0e-2\n" + 
                           " sfincs_Er_option='zero'\n" + 
                           "!  sfincs_Er_option='estimate'\n\n")
        self.OPTIMUM_PARAMS['Profile_Functions'].set(extra_lines_str)

        # do stuff and things
        # Make two frames. One for run parameters and one for scana
        self.optimum_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="OPTIMUM Parameters",
                                   font=('Helvetica', '14'))
        self.optimum_frame.grid(row=1,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        counter = 0
        for this_key in ('NFUNC_MAX', 'EQUIL_TYPE', 'BOOTCALC_TYPE',
                         'OPT_TYPE', 'FTOL', 'XTOL', 'GTOL', 'FACTOR',
                         'EPSFCN', 'MODE', 'LKEEP_MINS'):
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.OPTIMUM_PARAMS,
                                              self.optimum_frame)

        counter += 1
        self.OPTIMUM_EXTRA_LINES = self.add_text_entry_to_frame(counter,
                                     'Extra_Lines',
                                     self.OPTIMUM_PARAMS,
                                     self.optimum_frame, entry_width=132,
                                     entry_height = 9)



    def make_targets_tab(self, this_tab):
        # Store everything in a dictionary for later
        self.OPTIMUM_TARGETS_PARAMS = {}

        # regcoil settings
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_CHI2_B'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_CHI2_B']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_CHI2_B']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_CHI2_B']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_CHI2_B']['TARGET'].set(0.0)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_CHI2_B']['SIGMA'].set(1.0e-3)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_CHI2_B']['COUNT'].set(1)
        
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_BNORMAL_TOTAL'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_BNORMAL_TOTAL']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_BNORMAL_TOTAL']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_BNORMAL_TOTAL']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_BNORMAL_TOTAL']['TARGET'].set(0.0)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_BNORMAL_TOTAL']['SIGMA'].set(1.0e-3)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_BNORMAL_TOTAL']['COUNT'].set(4096)
         
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_BNORMAL'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_BNORMAL']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_BNORMAL']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_BNORMAL']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_BNORMAL']['TARGET'].set(0.0)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_BNORMAL']['SIGMA'].set(1.0e-3)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_BNORMAL']['COUNT'].set(1)
        
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K']['TARGET'].set(4.0e6)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K']['SIGMA'].set(1.0e6)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K']['COUNT'].set(1)

        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K']['TARGET'].set(2.4e6)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K']['SIGMA'].set(1.0e6)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K']['COUNT'].set(1)

        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_C2P_DIST_MIN'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_C2P_DIST_MIN']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_C2P_DIST_MIN']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_C2P_DIST_MIN']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_C2P_DIST_MIN']['TARGET'].set(0.20)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_C2P_DIST_MIN']['SIGMA'].set(1.0e-5)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_C2P_DIST_MIN']['COUNT'].set(1)

        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_VOLUME_COIL'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_VOLUME_COIL']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_VOLUME_COIL']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_VOLUME_COIL']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_VOLUME_COIL']['TARGET'].set(20.0)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_VOLUME_COIL']['SIGMA'].set(0.5)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_VOLUME_COIL']['COUNT'].set(1)

        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_LAMBDA'] = {}
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_LAMBDA']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_LAMBDA']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_LAMBDA']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_LAMBDA']['TARGET'].set(0.0)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_LAMBDA']['SIGMA'].set(1.0e9)
        self.OPTIMUM_TARGETS_PARAMS['REGCOIL_LAMBDA']['COUNT'].set(1)

        # cobravmec settings
        self.OPTIMUM_TARGETS_PARAMS['BALLOON'] = {}
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['BALLOON_THETA'] = tk.StringVar()
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['BALLOON_THETA'].set('0.0 45.0 90.0')
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['BALLOON_ZETA'] = tk.StringVar()
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['BALLOON_ZETA'].set('0.0 45.0 90.0')
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['TARGET'].set(0.0)
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['SIGMA'].set(1.0)
        self.OPTIMUM_TARGETS_PARAMS['BALLOON']['COUNT'].set(51)
        
        # Boozer Coordinates
        self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD']  = {}
        self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD']['MBOZ']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD']['MBOZ'].set(64)
        self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD']['NBOZ']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD']['NBOZ'].set(32)
        
        # Boozer Coordinate Helicity
        self.OPTIMUM_TARGETS_PARAMS['HELICITY'] = tk.StringVar()
        self.OPTIMUM_TARGETS_PARAMS['HELICITY'].set('\n')
        
        # NEO
        self.OPTIMUM_TARGETS_PARAMS['NEO'] = tk.StringVar()
        self.OPTIMUM_TARGETS_PARAMS['NEO'].set('\n')

        # MAGWELL
        self.OPTIMUM_TARGETS_PARAMS['MAGWELL'] = {}
        self.OPTIMUM_TARGETS_PARAMS['MAGWELL']['TARGET']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['MAGWELL']['SIGMA']  = tk.DoubleVar()
        self.OPTIMUM_TARGETS_PARAMS['MAGWELL']['COUNT']  = tk.IntVar()
        self.OPTIMUM_TARGETS_PARAMS['MAGWELL']['TARGET'].set(0.01)
        self.OPTIMUM_TARGETS_PARAMS['MAGWELL']['SIGMA'].set(-1.0e30)
        self.OPTIMUM_TARGETS_PARAMS['MAGWELL']['COUNT'].set(51)

        # GAMMA_C
        self.OPTIMUM_TARGETS_PARAMS['GAMMA_C'] = tk.StringVar()
        self.OPTIMUM_TARGETS_PARAMS['GAMMA_C'].set('\n')

        # Aspect Ratio
        self.OPTIMUM_TARGETS_PARAMS['ASPECT'] = tk.StringVar()
        self.OPTIMUM_TARGETS_PARAMS['ASPECT'].set('\n')

        # 'Optimum Extras'
        for this_key in ('BALLOON', 'BOOZER_COORD'):
            self.OPTIMUM_TARGETS_PARAMS[this_key]['Enabled'] = tk.BooleanVar()
            self.OPTIMUM_TARGETS_PARAMS[this_key]['Enabled'].set(False)


        # do stuff and things
        # Make frames for each target set. One for  run parameters and one for extras
        counter_frame = 0
        # regcoil
        counter_frame +=1
        self.optimum_regcoil_targets_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="REGCOIL (Regularized NESCoil)",
                                   font=('Helvetica', '14'))
        self.optimum_regcoil_targets_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        counter = 1
        self.add_optimum_legend_to_frame(counter, self.optimum_regcoil_targets_frame)
        
        for this_key in ('REGCOIL_CHI2_B', 'REGCOIL_BNORMAL_TOTAL', 'REGCOIL_MAX_BNORMAL', 'REGCOIL_MAX_K', 'REGCOIL_RMS_K',
                         'REGCOIL_C2P_DIST_MIN', 'REGCOIL_VOLUME_COIL', 'REGCOIL_LAMBDA'):
            counter += 1
            self.add_optimum_target_entry_to_frame(counter, this_key,
                                              self.OPTIMUM_TARGETS_PARAMS,
                                              self.optimum_regcoil_targets_frame)

        # cobravmec
        counter_frame +=1
        self.optimum_cobravmec_targets_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="COBRAVMEC (Ideal MHD Ballooning Stability)",
                                   font=('Helvetica', '14'))
        self.optimum_cobravmec_targets_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        counter = 1
        for this_key in ('BALLOON_THETA', 'BALLOON_ZETA'):
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.OPTIMUM_TARGETS_PARAMS['BALLOON'],
                                              self.optimum_cobravmec_targets_frame)

        counter += 1
        self.add_optimum_legend_to_frame(counter, self.optimum_cobravmec_targets_frame)

        counter += 1
        self.add_optimum_target_entry_to_frame(counter, 'BALLOON',
                                          self.OPTIMUM_TARGETS_PARAMS,
                                          self.optimum_cobravmec_targets_frame)

        
        # Boozer coordinates
        counter_frame +=1
        self.optimum_boozer_coordinates_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="BOOZ_XFORM (Boozer Coordinate Transformation)",
                                   font=('Helvetica', '14'))
        self.optimum_boozer_coordinates_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        counter = 1
        for this_key in ('MBOZ', 'NBOZ'):
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD'],
                                              self.optimum_boozer_coordinates_frame)


        # MAGWELL 
        counter_frame +=1
        self.optimum_magwell_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="MAGWELL",
                                   font=('Helvetica', '14'))
        self.optimum_magwell_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        counter = 1
        self.add_optimum_target_entry_to_frame(counter, 'MAGWELL',
                                          self.OPTIMUM_TARGETS_PARAMS,
                                          self.optimum_magwell_frame)


        # Helicity coordinates
        counter_frame = 1
        self.optimum_helicity_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="Helicity",
                                   font=('Helvetica', '14'))
        self.optimum_helicity_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        sfincs_file = open('templates/helicity_init_v0.txt')
        sfincs_text = sfincs_file.read()
        sfincs_file.close()
        self.OPTIMUM_TARGETS_PARAMS['HELICITY'].set(sfincs_text)

        counter +=1 # frame 3 will be used
        self.Helicity_Text = self.add_text_entry_to_frame(counter,
                                    'HELICITY',
                                    self.OPTIMUM_TARGETS_PARAMS,
                                    self.optimum_helicity_frame, 
                                    entry_width=80, entry_height=8)


        # NEO 
        counter_frame +=1
        self.optimum_neo_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="NEO",
                                   font=('Helvetica', '14'))
        self.optimum_neo_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        sfincs_file = open('templates/neo_init_v0.txt')
        sfincs_text = sfincs_file.read()
        sfincs_file.close()
        self.OPTIMUM_TARGETS_PARAMS['NEO'].set(sfincs_text)

        counter +=1 # 
        self.NEO_Text = self.add_text_entry_to_frame(counter,
                                    'NEO',
                                    self.OPTIMUM_TARGETS_PARAMS,
                                    self.optimum_neo_frame, 
                                    entry_width=80, entry_height=8)



        # Gamma_C 
        counter_frame +=1
        self.gammac_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="GAMMA_C",
                                   font=('Helvetica', '14'))
        self.gammac_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        sfincs_file = open('templates/gamma_c_init_v0.txt')
        sfincs_text = sfincs_file.read()
        sfincs_file.close()
        self.OPTIMUM_TARGETS_PARAMS['GAMMA_C'].set(sfincs_text)

        counter +=1 # 
        self.GAMMA_C_Text = self.add_text_entry_to_frame(counter,
                                    'GAMMA_C',
                                    self.OPTIMUM_TARGETS_PARAMS,
                                    self.gammac_frame, 
                                    entry_width=80, entry_height=5)


        # Aspect Ratio 
        counter_frame +=1
        self.optimum_aspect_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="Aspect Ratio",
                                   font=('Helvetica', '14'))
        self.optimum_aspect_frame.grid(row=counter_frame,
                        rowspan=1,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        sfincs_file = open('templates/aspect_init_v0.txt')
        sfincs_text = sfincs_file.read()
        sfincs_file.close()
        self.OPTIMUM_TARGETS_PARAMS['ASPECT'].set(sfincs_text)

        counter +=1 # 
        self.ASPECT_Text = self.add_text_entry_to_frame(counter,
                                    'ASPECT',
                                    self.OPTIMUM_TARGETS_PARAMS,
                                    self.optimum_aspect_frame, 
                                    entry_width=80, entry_height=4)












    def make_variables_tab(self, this_tab):
        # Store everything in a dictionary for later
        self.OPTIMUM_VARIABLES_PARAMS = {}
        # These parameters need to be written out later.
        self.OPTIMUM_VARIABLES_PARAMS['Variable_Lines'] = tk.StringVar()
        extra_lines_str = ("")
        
        # type of tk.Text
        self.OPTIMUM_VARIABLES_PARAMS['Variable_Lines'].set(extra_lines_str)

        # do stuff and things
        # Make two frames. One for parameters and one for 'append buttons'
        self.optimum_variables_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="VARIABLES",
                                   font=('Helvetica', '14'))
        self.optimum_variables_frame.grid(row=1,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        self.optimum_append_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="Append Options",
                                   font=('Helvetica', '14'))
        self.optimum_append_frame.grid(row=2,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        counter = 0
        for this_key in ():
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.OPTIMUM_VARIABLES_PARAMS,
                                              self.optimum_variables_frame)

         
        self.OPTIMUM_VARIABLE_LINES = self.add_text_entry_to_frame(counter,
                                     'Variable_Lines',
                                     self.OPTIMUM_VARIABLES_PARAMS,
                                     self.optimum_variables_frame, entry_width=132,
                                     entry_height = 30)

        row_counter = 1
        col_counter = 1
        
        # REGCOIL OPT
        clear_variables_button = tk.Button(self.optimum_append_frame,
                                command=self.clear_variables,
                                text='Clear Variables')
        clear_variables_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        row_counter = 2
        col_counter = 1
        regcoil_0_button = tk.Button(self.optimum_append_frame,
                                command=None, # self.doit,
                                text='Regcoil, Uniform')
        regcoil_0_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_1_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_2x2, # self.doit,
                                text='Regcoil, 2x2')
        regcoil_1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_2_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_4x4, # self.doit,
                                text='Regcoil, 4x4')
        regcoil_2_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_3_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_6x6, # self.doit,
                                text='Regcoil, 6x6')
        regcoil_3_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_4_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_8x8, # self.doit,
                                text='Regcoil, 8x8')
        regcoil_4_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_5_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_10x10, # self.doit,
                                text='Regcoil, 10x10')
        regcoil_5_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_6_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_12x12, # self.doit,
                                text='Regcoil, 12x12')
        regcoil_6_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_7_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_16x16, # self.doit,
                                text='Regcoil, 16x16')
        regcoil_7_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        regcoil_8_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_regcoil_24x16, # self.doit,
                                text='Regcoil, 24x16')
        regcoil_8_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)


        # Boundary opt
        row_counter = 3
        col_counter = 1
        boundary_1_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_boundary_2x2, # self.doit,
                                text='Boundary, 2x2')
        boundary_1_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        boundary_2_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_boundary_4x4, # self.doit,
                                text='Boundary, 4x4')
        boundary_2_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        boundary_3_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_boundary_6x6, # self.doit,
                                text='Boundary, 6x6')
        boundary_3_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)

        col_counter += 1
        boundary_4_button = tk.Button(self.optimum_append_frame,
                                command=self.append_variables_boundary_8x8, # self.doit,
                                text='Boundary, 8x8')
        boundary_4_button.grid(row=row_counter, rowspan=1, column=col_counter, columnspan=1)




    def clear_variables(self):
        my_empty_tk_string = tk.StringVar()
        my_empty_tk_string.set('')
        self.OPTIMUM_VARIABLE_LINES.delete('1.0', 'end')
        self.OPTIMUM_VARIABLE_LINES.insert('1.0', my_empty_tk_string.get())


    # REGCOIL Files
    def append_variables_regcoil_uniform(self):
        file_handle = open('templates/append_regcoil_uniform.txt');
        file_text = read(file_handle)
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.delete('1.0', 'end')
        self.OPTIMUM_VARIABLE_LINES.insert('1.0', new_text.get())

    def append_variables_regcoil_2x2(self):
        file_handle = open('templates/append_regcoil_2x2.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_regcoil_4x4(self):
        file_handle = open('templates/append_regcoil_4x4.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_regcoil_6x6(self):
        file_handle = open('templates/append_regcoil_6x6.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_regcoil_8x8(self):
        file_handle = open('templates/append_regcoil_8x8.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_regcoil_10x10(self):
        file_handle = open('templates/append_regcoil_10x10.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_regcoil_12x12(self):
        file_handle = open('templates/append_regcoil_12x12.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_regcoil_16x16(self):
        file_handle = open('templates/append_regcoil_16x16.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_regcoil_24x16(self):
        file_handle = open('templates/append_regcoil_24x16.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())


    # Boundary files
    def append_variables_boundary_2x2(self):
        file_handle = open('templates/append_boundary_2x2.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_boundary_4x4(self):
        file_handle = open('templates/append_boundary_4x4.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_boundary_6x6(self):
        file_handle = open('templates/append_boundary_6x6.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

    def append_variables_boundary_8x8(self):
        file_handle = open('templates/append_boundary_8x8.txt');
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())


    def make_regcoil_tab(self, this_tab):
        # VMEC Execution parameters.
        # Store everything in a dictionary for later
        self.REGCOIL_PARAMS = {}

        self.REGCOIL_PARAMS['nlambda'] = tk.IntVar()
        self.REGCOIL_PARAMS['nlambda'].set(51)

        self.REGCOIL_PARAMS['ntheta_plasma'] = tk.IntVar()
        self.REGCOIL_PARAMS['ntheta_plasma'].set(64)

        self.REGCOIL_PARAMS['ntheta_coil'] = tk.IntVar()
        self.REGCOIL_PARAMS['ntheta_coil'].set(64)

        self.REGCOIL_PARAMS['nzeta_plasma'] = tk.IntVar()
        self.REGCOIL_PARAMS['nzeta_plasma'].set(64)

        self.REGCOIL_PARAMS['nzeta_coil'] = tk.IntVar()
        self.REGCOIL_PARAMS['nzeta_coil'].set(64)

        self.REGCOIL_PARAMS['mpol_potential'] = tk.IntVar()
        self.REGCOIL_PARAMS['mpol_potential'].set(12)

        self.REGCOIL_PARAMS['ntor_potential'] = tk.IntVar()
        self.REGCOIL_PARAMS['ntor_potential'].set(12)

        self.REGCOIL_PARAMS['general_option'] = tk.IntVar()
        self.REGCOIL_PARAMS['general_option'].set(5)

        self.REGCOIL_PARAMS['geometry_option_plasma'] = tk.IntVar()
        self.REGCOIL_PARAMS['geometry_option_plasma'].set(2)

        self.REGCOIL_PARAMS['geometry_option_coil'] = tk.IntVar()
        self.REGCOIL_PARAMS['geometry_option_coil'].set(3)


        self.REGCOIL_PARAMS['net_poloidal_current_Amperes'] = tk.DoubleVar()
        self.REGCOIL_PARAMS['net_poloidal_current_Amperes'].set(0.0)

        self.REGCOIL_PARAMS['net_toroidal_current_Amperes'] = tk.DoubleVar()
        self.REGCOIL_PARAMS['net_toroidal_current_Amperes'].set(0.0)

        self.REGCOIL_PARAMS['symmetry_option'] = tk.IntVar()
        self.REGCOIL_PARAMS['symmetry_option'].set(1)

        self.REGCOIL_PARAMS['target_option'] = tk.StringVar()
        self.REGCOIL_PARAMS['target_option'].set('max_K')


        self.REGCOIL_PARAMS['nescin_filename'] = tk.StringVar()
        self.REGCOIL_PARAMS['nescin_filename'].set('nescin.out')

        self.REGCOIL_PARAMS['load_bnorm'] = tk.BooleanVar()
        self.REGCOIL_PARAMS['load_bnorm'].set(tk.FALSE)

        self.REGCOIL_PARAMS['bnorm_filename'] = tk.StringVar()
        self.REGCOIL_PARAMS['bnorm_filename'].set('bnorm.in')


        # do stuff and things
        # Make two frames. One for VMEC run parameters and one for V3FIT
        self.regcoil_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="REGCOIL Parameters",
                                   font=('Helvetica', '14'))
        self.regcoil_frame.grid(row=1,
                        rowspan=2,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        # Loop over all of the VMEC items and add them to the VMEC frame
        # vmec_keys = self.REGCOIL_PARAMS.keys()
        # Not looping because the order changes at time. Listing the explicetly
        counter = 0
        for this_key in ('nlambda', 'ntheta_plasma', 'ntheta_coil', 'nzeta_plasma', 'nzeta_coil',
                         'mpol_potential', 'ntor_potential', 'general_option', 'geometry_option_plasma', 
                         'geometry_option_coil',
                         'target_option', 'nescin_filename', 'net_poloidal_current_Amperes', 
                         'net_toroidal_current_Amperes', 'symmetry_option', 'load_bnorm',
                         'bnorm_filename'):
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.REGCOIL_PARAMS,
                                              self.regcoil_frame)
          
    def make_bootstrap_tab(self, this_tab):
        # VMEC Execution parameters.
        # Store everything in a dictionary for later
        # BOOTSJ
        self.BOOTSTRAP_PARAMS = {}
        self.BOOTSTRAP_PARAMS['BOOTSJ'] = {}
        self.BOOTSTRAP_PARAMS['BOOTSJ']['Enabled'] = tk.BooleanVar()
        self.BOOTSTRAP_PARAMS['BOOTSJ']['Enabled'].set(False)
        self.BOOTSTRAP_PARAMS['MBUSE']  = tk.IntVar()
        self.BOOTSTRAP_PARAMS['MBUSE'].set(64)
        self.BOOTSTRAP_PARAMS['NBUSE']  = tk.IntVar()
        self.BOOTSTRAP_PARAMS['NBUSE'].set(32)
        self.BOOTSTRAP_PARAMS['ZEFF1']  = tk.DoubleVar()
        self.BOOTSTRAP_PARAMS['ZEFF1'].set(1.00)
        self.BOOTSTRAP_PARAMS['DENS0']  = tk.DoubleVar()
        self.BOOTSTRAP_PARAMS['DENS0'].set(0.70)
        self.BOOTSTRAP_PARAMS['TETI']  = tk.DoubleVar()
        self.BOOTSTRAP_PARAMS['TETI'].set(1.00)
        self.BOOTSTRAP_PARAMS['TEMPRES']  = tk.DoubleVar()
        self.BOOTSTRAP_PARAMS['TEMPRES'].set(-1.0)
        self.BOOTSTRAP_PARAMS['DAMP_BS']  = tk.DoubleVar()
        self.BOOTSTRAP_PARAMS['DAMP_BS'].set(0.1)
        self.BOOTSTRAP_PARAMS['ISYMM0']  = tk.IntVar()
        self.BOOTSTRAP_PARAMS['ISYMM0'].set(0)
        self.BOOTSTRAP_PARAMS['ATE']  = tk.StringVar()
        self.BOOTSTRAP_PARAMS['ATE'].set('2.0 -2.0')
        self.BOOTSTRAP_PARAMS['ATI']  = tk.StringVar()
        self.BOOTSTRAP_PARAMS['ATI'].set('2.0 -2.0')

        # SFINCS Stuff
        self.BOOTSTRAP_PARAMS['SFINCS'] = {}
        self.BOOTSTRAP_PARAMS['SFINCS']['Enabled'] = tk.BooleanVar()
        self.BOOTSTRAP_PARAMS['SFINCS']['Enabled'].set(False)
        self.BOOTSTRAP_PARAMS['SFINCS_Text'] = tk.StringVar()
        sfincs_file = open('templates/sfincs_init_v0.txt')
        sfincs_text = sfincs_file.read()
        sfincs_file.close()
        self.BOOTSTRAP_PARAMS['SFINCS_Text'].set(sfincs_text)

        # Profiles Stuff
        self.BOOTSTRAP_PARAMS['PROFILES'] = {}
        self.BOOTSTRAP_PARAMS['PROFILES']['Enabled'] = tk.BooleanVar()
        self.BOOTSTRAP_PARAMS['PROFILES']['Enabled'].set(False)
        self.BOOTSTRAP_PARAMS['PROFILES_Text'] = tk.StringVar()
        sfincs_file = open('templates/profiles_init_v0.txt')
        sfincs_text = sfincs_file.read()
        sfincs_file.close()
        self.BOOTSTRAP_PARAMS['PROFILES_Text'].set(sfincs_text)


        # do stuff and things
        # Make two frames. 

        profiles_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="Profiles",
                                   font=('Helvetica', '14'))
        profiles_frame.grid(row=1,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        bootsj_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="BOOTSJ Parameters",
                                   font=('Helvetica', '14'))
        bootsj_frame.grid(row=1,
                        rowspan=1,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        sfincs_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="SFINCS Parameters",
                                   font=('Helvetica', '14'))
        sfincs_frame.grid(row=3,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)


        counter = 0
        counter2 = 1
        self.add_selectables_line_to_frame(counter, counter2, 'BOOTSJ',
                                           self.BOOTSTRAP_PARAMS, bootsj_frame)

        for this_key in ('MBUSE', 'NBUSE', 'ZEFF1', 'DENS0', 'TETI', 'TEMPRES',
                         'DAMP_BS', 'ISYMM0', 'ATE', 'ATI'):
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.BOOTSTRAP_PARAMS,
                                              bootsj_frame)

        row_counter = 0 # frame 2 will be used
        col_counter = 0
        self.add_selectables_line_to_frame(row_counter, col_counter, 'SFINCS',
                                           self.BOOTSTRAP_PARAMS, sfincs_frame)
        row_counter += 1 # frame 2 will be used
        self.SFINCS_Text = self.add_text_entry_to_frame(row_counter,
                                    'SFINCS_Text',
                                    self.BOOTSTRAP_PARAMS,
                                    sfincs_frame, entry_width=80)


        row_counter = 0 # frame 2 will be used
        col_counter = 0
        self.add_selectables_line_to_frame(row_counter, col_counter, 'PROFILES',
                                           self.BOOTSTRAP_PARAMS, profiles_frame)
        row_counter += 1 # frame 3 will be used
        self.Profiles_Text = self.add_text_entry_to_frame(row_counter,
                                    'PROFILES_Text',
                                    self.BOOTSTRAP_PARAMS,
                                    profiles_frame, entry_width=80)


        

    def make_scanables_tab(self, this_tab):
        # Store everything in a dictionary for later
        self.OPTIMUM_SCAN_PARAMS = {}
        #  Options with sigmas, or weights that are synced individually
        for this_key in ('REGCOIL_MAX_K', 'REGCOIL_RMS_K', 'REGCOIL_CHI2_B', 
                         'REGCOIL_BNORMAL_TOTAL', 'REGCOIL_MAX_BNORMAL', 'REGCOIL_VOLUME_COIL',
                        'REGCOIL_C2P_DIST_MIN', 'REGCOIL_LAMBDA'):
            self.OPTIMUM_SCAN_PARAMS[this_key] = {}
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'] = tk.BooleanVar()
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].set(False)
            self.OPTIMUM_SCAN_PARAMS[this_key]['Values'] = tk.StringVar()
            self.OPTIMUM_SCAN_PARAMS[this_key]['Values'].set('')
            self.OPTIMUM_SCAN_PARAMS[this_key]['Sigmas'] = tk.StringVar()
            self.OPTIMUM_SCAN_PARAMS[this_key]['Sigmas'].set('')
            
        # Options that don't have associated weights
        for this_key in ('FACTOR', 'EPSFCN'):
            self.OPTIMUM_SCAN_PARAMS[this_key] = {}
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'] = tk.BooleanVar()
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].set(False)
            self.OPTIMUM_SCAN_PARAMS[this_key]['Values'] = tk.StringVar()
            self.OPTIMUM_SCAN_PARAMS[this_key]['Values'].set('')

        # 'Regcoil Fourier Selectables'
        for this_key in ('REGCOIL_FOURIER_SPECTRUM', 'REGCOIL_2x2', 'REGCOIL_4x4', 'REGCOIL_6x6', 'REGCOIL_8x8',
                         'REGCOIL_10x10', 'REGCOIL_12x12', 'REGCOIL_16x16', 'REGCOIL_24x16',
                         'COORDINATE_REGCOIL_TARGET_VALUE', 'AUTOGEN_REGCOIL_BOUNDS', 'AUTOGEN_REGCOIL_D',
                         'CLEAR_VARS_B4_REGCOIL', 'CLEAR_VARS_B4_BOUNDARY'):
            self.OPTIMUM_SCAN_PARAMS[this_key] = {}
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'] = tk.BooleanVar()
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].set(False)

        # 'Boundary Fourier Selectables'
        for this_key in ('BOUNDARY_SPECTRUM', 'BOUND_2x2', 'BOUND_4x4', 'BOUND_6x6', 'BOUND_8x8',
                         'AUTOGEN_BOUNDARY_LIMITS', 'AUTOGEN_BOUNDARY_D'):
            self.OPTIMUM_SCAN_PARAMS[this_key] = {}
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'] = tk.BooleanVar()
            self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].set(False)

        self.optimum_scan_frame = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="Scan-able Parameters",
                                   font=('Helvetica', '14'))

        self.optimum_scan_frame.grid(row=2,
                        rowspan=1,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)


        # make max k / rms k selection - left to the user
        # make lm / lbfgsb selection - left to the user
        doit_button = tk.Button(self.optimum_scan_frame,
                                command=self.doit, # self.doit,
                                text='<==-- Do it --==>')

        doit_button.grid(row=1, rowspan=1, column=1, columnspan=1)
        
        counter = 1  # row counter          
        counter2 = 0 # column counter

        this_key = 'COORDINATE_REGCOIL_TARGET_VALUE'
        counter += 1
        counter2 += 1
        self.add_selectables_line_to_frame(counter, counter2, this_key, self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)
        
        for this_key in ('REGCOIL_MAX_K', 'REGCOIL_RMS_K', 'REGCOIL_CHI2_B', 'REGCOIL_BNORMAL_TOTAL', 'REGCOIL_MAX_BNORMAL',
                         'REGCOIL_VOLUME_COIL', 'REGCOIL_C2P_DIST_MIN', 'REGCOIL_LAMBDA'):
            counter += 1
            self.add_label_and_scan_line_to_frame(counter, this_key,
                                              self.OPTIMUM_SCAN_PARAMS,
                                              self.optimum_scan_frame,
                                              with_weights=True, entry_width=45, col_span=3)
        for this_key in ('FACTOR', 'EPSFCN'):
            counter += 1
            self.add_label_and_scan_line_to_frame(counter, this_key,
                                              self.OPTIMUM_SCAN_PARAMS,
                                              self.optimum_scan_frame,
                                              with_weights=False, entry_width=45, col_span=3)
                 
        counter += 1
        counter2 = 0    # column counter  
        for this_key in ('REGCOIL_FOURIER_SPECTRUM', 'REGCOIL_2x2', 'REGCOIL_4x4', 'REGCOIL_6x6', 'REGCOIL_8x8',
                         'REGCOIL_10x10', 'REGCOIL_12x12', 'REGCOIL_16x16', 'REGCOIL_24x16'):
            counter2 += 1
            self.add_selectables_line_to_frame(counter, counter2, this_key, self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)
              
        counter +=1
        counter2 = 0
        counter2 += 1
        self.add_selectables_line_to_frame(counter, counter2, 'CLEAR_VARS_B4_REGCOIL', self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)

        counter += 0
        counter2 += 1
        self.add_selectables_line_to_frame(counter, counter2, 'CLEAR_VARS_B4_BOUNDARY', self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)

        counter +=1
        counter2 = 0
        counter2 += 1
        self.add_selectables_line_to_frame(counter, counter2, 'AUTOGEN_REGCOIL_BOUNDS', self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)

        counter += 0
        counter2 += 1
        self.add_selectables_line_to_frame(counter, counter2, 'AUTOGEN_REGCOIL_D', self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)
        
        counter += 1
        counter2 = 0    # column counter  
        for this_key in ('BOUNDARY_SPECTRUM', 'BOUND_2x2', 'BOUND_4x4', 'BOUND_6x6', 'BOUND_8x8'):
            counter2 += 1
            self.add_selectables_line_to_frame(counter, counter2, this_key, self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)
              
        counter +=1
        counter2 = 0
        counter2 += 1
        self.add_selectables_line_to_frame(counter, counter2, 'AUTOGEN_BOUNDARY_LIMITS', self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)

        counter2 += 1
        self.add_selectables_line_to_frame(counter, counter2, 'AUTOGEN_BOUNDARY_D', self.OPTIMUM_SCAN_PARAMS, self.optimum_scan_frame)


    def make_filenames_tab(self, this_tab):
        # Store everything in a dictionary for later
        self.FILESETC = {} # FILESETC = "Files, etc"

        self.FILESETC['file_base'] = tk.StringVar()
        self.FILESETC['file_base'].set('input.stell0')

        self.FILESETC['folder_base'] = tk.StringVar()
        self.FILESETC['folder_base'].set('stell0')

        self.FILESETC['README'] = tk.StringVar()
        extra_lines_str = (' ! Contents of the README file go here\n' +
                           ' ! \n')
        self.FILESETC['README'].set(extra_lines_str)

        self.FILESETC['chtc'] = tk.StringVar()
        self.FILESETC['chtc'].set('chtc-stellopt.sub')

        self.FILESETC['chtc_contents'] = tk.StringVar()
        extra_lines_str = ("# stellopt-cthc.sub\n" +
                           "# starter submit file for CHTC jobs\n" +
                           " \n" +
                           " universe = docker\n" +
                           " docker_image = benjaminfaber/stellopt_wisc:regcoil_2019a\n" +
                           " \n"                       +
                           " log = job_$(Cluster).log\n" +
                           " error = job_$(Cluster)_$(Process).err\n" +
                           " output = job_$(Cluster)_$(Process).out\n" +
                           " \n" +
                           " executable = submit_stellopt.sh\n"  +
                           " arguments =\n" +
                           " requirements = (TARGET.has_avx == True)\n"  +
                           " \n" +
                           " should_transfer_files = YES\n" +
                           " when_to_transfer_output = ON_EXIT\n" +
                           " transfer_input_files = input.stell0,nescin.out\n" +
                           " \n"                                +
                           " \n" +
                           " request_cpus = 8\n" +
                           " request_memory = 16GB\n" +
                           " \n" +
                           " request_disk =  10GB\n"  +
                           " \n" +
                           " queue 1 \n")
        self.FILESETC['chtc_contents'].set(extra_lines_str)

        self.FILESETC['submit'] = tk.StringVar()
        self.FILESETC['submit'].set('submit_stellopt.sh')

        self.FILESETC['submit_contents'] = tk.StringVar()
        extra_lines_str = ('#!/bin/sh\n\n' +
                           'mpiexec -n 8 --allow-run-as-root ' +
                           '/bin/xstelloptv2 input.stell0\n' + 
                           'rm wout_stell0_opt* regcoil_nescout.stell0_opt*')
        self.FILESETC['submit_contents'].set(extra_lines_str)

        self.FILESETC['NESCIN'] = {}
        self.FILESETC['NESCIN']['Enabled'] = tk.BooleanVar()
        self.FILESETC['NESCIN']['Enabled'].set(False)
        self.FILESETC['NESCIN_FILEIN'] = tk.StringVar()
        self.FILESETC['NESCIN_FILEIN'].set('')
        self.FILESETC['NESCIN_FILEOUT'] = tk.StringVar()
        self.FILESETC['NESCIN_FILEOUT'].set('nescin.out')
        
        self.FILESETC['BOUNDARY_INIT'] = {}
        self.FILESETC['BOUNDARY_INIT']['Enabled'] = tk.BooleanVar()
        self.FILESETC['BOUNDARY_INIT']['Enabled'].set(False)
        self.FILESETC['BOUNDARY_INIT_IN'] = tk.StringVar()
        self.FILESETC['BOUNDARY_INIT_IN'].set('/Users/schmittj/src/xStellTools/Generator/templates/vmec_mljs4_a3b25_p2_nml')

        
        # Make two frames. 
        self.filenames_frame1 = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="Folders & Filenames",
                                   font=('Helvetica', '14'))
        self.filenames_frame1.grid(row=1,
                        rowspan=2,
                        column=1,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        self.filenames_frame2 = tk.LabelFrame(this_tab,
                                   bg=self.bg_color_1,
                                   bd=5,
                                   padx=5,
                                   pady=5,
                                   relief=tk.RIDGE,
                                   text="Additional Filenames",
                                   font=('Helvetica', '14'))
        self.filenames_frame2.grid(row=1,
                        rowspan=2,
                        column=2,
                        columnspan=1,
                        padx=5, pady=5, ipadx=5, ipady=5)

        # Loop over all of the VMEC items and add them to the VMEC frame
        # vmec_keys = self.REGCOIL_PARAMS.keys()
        # Not looping because the order changes at time. Listing the explicetly
        counter = 0
        for this_key in ('file_base', 'folder_base'):
            counter += 1
            self.add_label_and_entry_to_frame(counter, this_key,
                                              self.FILESETC,
                                              self.filenames_frame1)
        counter += 1
        self.README_LINES = self.add_text_entry_to_frame(counter,
                                     'README',
                                     self.FILESETC,
                                     self.filenames_frame1, entry_width=80,
                                     entry_height = 8)
        counter += 1
        self.add_label_and_entry_to_frame(counter, 'chtc',
                                          self.FILESETC,
                                          self.filenames_frame1)
        counter += 1
        self.CHTC_CONTENTS = self.add_text_entry_to_frame(counter,
                                     'chtc_contents',
                                     self.FILESETC,
                                     self.filenames_frame1, entry_width=80,
                                     entry_height = 28)
        counter += 1
        self.add_label_and_entry_to_frame(counter, 'submit',
                                          self.FILESETC,
                                          self.filenames_frame1)
        counter += 0
        select_stellopt_submit = tk.Button(self.filenames_frame1,
                                command=self.find_stellopt_submit, # self.doit,
                                text='Select stellopt submit file')

        select_stellopt_submit.grid(row=counter, rowspan=1, column=2, columnspan=1)
        
        counter += 1
        self.SUBMIT_CONTENTS = self.add_text_entry_to_frame(counter,
                                     'submit_contents',
                                     self.FILESETC,
                                     self.filenames_frame1, entry_width=80,
                                     entry_height = 6)
        
        counter = 0
        counter2 = 0
        
        counter2 = 1
        self.add_selectables_line_to_frame(counter, counter2, 'NESCIN',
                                           self.FILESETC, self.filenames_frame2)
        counter += 1
        self.add_label_and_entry_to_frame(counter, 'NESCIN_FILEIN',
                                          self.FILESETC,
                                          self.filenames_frame2)
        counter += 1
        self.add_label_and_entry_to_frame(counter, 'NESCIN_FILEOUT',
                                          self.FILESETC,
                                          self.filenames_frame2)

        counter += 1
        find_nescin = tk.Button(self.filenames_frame2,
                                command=self.find_nescin, # self.doit,
                                text='Select NESCIN File')

        find_nescin.grid(row=counter, rowspan=1, column=1, columnspan=1)


        
        counter2 = 1
        counter += 1
        self.add_selectables_line_to_frame(counter, counter2, 'BOUNDARY_INIT',
                                           self.FILESETC, self.filenames_frame2)
        counter += 1
        self.add_label_and_entry_to_frame(counter, 'BOUNDARY_INIT_IN',
                                          self.FILESETC,
                                          self.filenames_frame2)

        counter += 1
        find_boundary = tk.Button(self.filenames_frame2,
                                command=self.find_boundary, # self.doit,
                                text='Select BOUNDARY INIT File')

        find_boundary.grid(row=counter, rowspan=1, column=1, columnspan=1)
        
    def find_nescin(self):
        # Use input as the intial guess
        importFile = filedialog.askopenfilename(
            initialdir=os.path.dirname('./input/'))
        self.FILESETC['NESCIN_FILEIN'].set(importFile)

    def find_boundary(self):
        # Use input as the intial guess
        importFile = filedialog.askopenfilename(
            initialdir=os.path.dirname('./templates/'))
        self.FILESETC['BOUNDARY_INIT_IN'].set(importFile)
       
    def find_stellopt_submit(self):
        import_file = filedialog.askopenfilename(
            initialdir=os.path.dirname('./input/'))
        (the_folder, the_file) = os.path.split(import_file)
        self.FILESETC['submit'].set(the_file)
        file_handle = open(import_file);
        file_text = file_handle.read()
        file_handle.close()
        new_text = tk.StringVar()
        new_text.set(file_text)
        self.SUBMIT_CONTENTS.delete('1.0', 'end')
        self.SUBMIT_CONTENTS.insert('1.0', new_text.get())
          
    def doit(self):
        # first, count the # of variations and determine the total nubmer of files to be generated
        total_variations = 1
        self.variation_list = []
        self.variation_count = []
        for this_key in ('REGCOIL_MAX_K', 'REGCOIL_RMS_K', 'REGCOIL_CHI2_B', 'REGCOIL_BNORMAL_TOTAL', 'REGCOIL_MAX_BNORMAL',
                         'REGCOIL_VOLUME_COIL', 'REGCOIL_C2P_DIST_MIN', 'REGCOIL_LAMBDA'):
            if (self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].get()):
                selected_values = np.array(self.OPTIMUM_SCAN_PARAMS[this_key]['Values'].get().replace(',', ' ').split(), 'float')
                selected_sigmas = np.array(self.OPTIMUM_SCAN_PARAMS[this_key]['Sigmas'].get().replace(',', ' ').split(), 'float')
                print(str(selected_values))
                print(str(selected_sigmas))
                
                if ( len(selected_values) != len(selected_sigmas) ):
                    print('<----Different number of values and sigmas: ' + this_key)
                    variations = 0
                else:
                    variations = len(selected_values)
                    self.OPTIMUM_SCAN_PARAMS[this_key]['values'] = selected_values
                    self.OPTIMUM_SCAN_PARAMS[this_key]['sigmas'] = selected_sigmas
                    self.variation_list.append(this_key)
                    self.variation_count.append(variations)
                total_variations = total_variations * variations
            else:
                variations = 1
                # not used?
                self.OPTIMUM_SCAN_PARAMS[this_key]['values'] = self.OPTIMUM_TARGETS_PARAMS[this_key]['TARGET'].get()
                self.OPTIMUM_SCAN_PARAMS[this_key]['sigmas'] = self.OPTIMUM_TARGETS_PARAMS[this_key]['SIGMA'].get()
            self.OPTIMUM_SCAN_PARAMS[this_key]['count'] = self.OPTIMUM_TARGETS_PARAMS[this_key]['COUNT'].get()
            self.OPTIMUM_SCAN_PARAMS[this_key]['variations'] = variations
        for this_key in ('FACTOR', 'EPSFCN'):
            if (self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].get()):
                selected_values = np.array(self.OPTIMUM_SCAN_PARAMS[this_key]['Values'].get().replace(',', ' ').split(), 'float')
                variations = len(selected_values)
                total_variations = total_variations * variations
                self.OPTIMUM_SCAN_PARAMS[this_key]['values'] = selected_values
                self.variation_list.append(this_key)
                self.variation_count.append(variations)
            else:
                # not used?
                variations = 1
                self.OPTIMUM_SCAN_PARAMS[this_key]['values'] = self.OPTIMUM_PARAMS[this_key].get()
            self.OPTIMUM_SCAN_PARAMS[this_key]['variations'] = variations

        if (self.OPTIMUM_SCAN_PARAMS['REGCOIL_FOURIER_SPECTRUM']['Enabled'].get() is True):
            self.variation_list.append('REGCOIL_FOURIER_SPECTRUM')
            num_rfs_selections = 0
            self.OPTIMUM_SCAN_PARAMS['REGCOIL_FOURIER_SPECTRUM']['Selections'] = []
            for this_key in ('REGCOIL_2x2', 'REGCOIL_4x4', 'REGCOIL_6x6', 'REGCOIL_8x8',
                             'REGCOIL_10x10', 'REGCOIL_12x12', 'REGCOIL_16x16', 'REGCOIL_24x16'):
                if (self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].get() is True):
                    num_rfs_selections += 1
                    self.OPTIMUM_SCAN_PARAMS['REGCOIL_FOURIER_SPECTRUM']['Selections'].append(this_key)
            # not used?
            variations = num_rfs_selections
            total_variations = total_variations * variations
            self.variation_count.append(variations)
            self.OPTIMUM_SCAN_PARAMS['REGCOIL_FOURIER_SPECTRUM']['variations'] = variations

 
        if (self.OPTIMUM_SCAN_PARAMS['BOUNDARY_SPECTRUM']['Enabled'].get() is True):
            self.variation_list.append('BOUNDARY_SPECTRUM')
            num_lcfs_selections = 0
            self.OPTIMUM_SCAN_PARAMS['BOUNDARY_SPECTRUM']['Selections'] = []
            for this_key in ('BOUND_2x2', 'BOUND_4x4', 'BOUND_6x6', 'BOUND_8x8'):
                if (self.OPTIMUM_SCAN_PARAMS[this_key]['Enabled'].get() is True):
                    num_lcfs_selections += 1
                    self.OPTIMUM_SCAN_PARAMS['BOUNDARY_SPECTRUM']['Selections'].append(this_key)
            # not used?
            variations = num_lcfs_selections
            total_variations = total_variations * variations
            self.variation_count.append(variations)
            self.OPTIMUM_SCAN_PARAMS['BOUNDARY_SPECTRUM']['variations'] = variations


        # confirm that these files should be generated
        if (total_variations == 0):
            tk.messagebox.showerror("Error", "Something is amiss...")
        else:
            verify_generate = tk.messagebox.askyesno("Question",
                    "Are you certain you want to generate "
                    + str(total_variations) + " stellopt job(s)?")
            if (verify_generate):
                print('<----OK!')
                self.generate_files(variations)
            else:
                print('<----It is always good to double-check!')
                
        return 
    
    
    def generate_files(self, variations):
        # generate files
        # generate the namelists for VMEC, OPTIMUM, REGCOIL, + more
        # Scans are handled inside the section for each separat namelist
        
        # At some point, it would be nice to use the f90nml package exclusively.
        # For now, things are mixed.
        # 
        # self.the_nml will be constructed to contain all of the lines that will be written
        # self.the_foldername will be constructed to be unique to the scan
        # self.the_filename will also be constructed (if neccessary)
        # self.my_tempfile 
        
        # initialize variables
        self.nml = {}
        self.the_foldername = self.FILESETC['folder_base'].get()
        self.the_filename = self.FILESETC['file_base'].get()
        self.my_tempfile1 = 'temp_sg1.out'
        self.my_tempfile2 = 'temp_sg2.out'
        self.my_tempfile3 = 'temp_sg3.out'
        self.my_tempfile4 = 'temp_sg4.out'
        self.my_tempfile5 = 'temp_sg5.out'

        self.recursive_scan_variation(self.variation_list,
                                      self.variation_count, 
                                      self.the_foldername, 
                                      self.the_filename)
        return

            
    def recursive_scan_variation(self, variation_list_in, 
                                 variation_count_in, 
                                 the_foldername_in,
                                 the_filename_in):
        # if no items to pop, make the file!
        if (len(variation_list_in) == 0):
          #        check self-consistent req's
          #        check self-consistent req's
          #        check self-consistent req's
          #        check self-consistent req's
          #        write namelists
          #        create files
          #        
            self.write_indata_nml()
            self.write_optimum_nml()
            self.write_regcoil_nml()
            if (self.BOOTSTRAP_PARAMS['BOOTSJ']['Enabled'] is True):
               self.write_bootsj_nml()
            if (self.BOOTSTRAP_PARAMS['SFINCS']['Enabled'] is True):
               self.write_sfincs_nml()
            self.create_files(the_foldername_in,the_filename_in)
        else: 
            # pop a variation of the list
            my_scan_item = variation_list_in.pop()
            my_scan_count = variation_count_in.pop()
            # depending on the type, do different things
            if (my_scan_item in('REGCOIL_MAX_K', 'REGCOIL_RMS_K', 'REGCOIL_CHI2_B', 'REGCOIL_BNORMAL_TOTAL', 'REGCOIL_MAX_BNORMAL',
                         'REGCOIL_VOLUME_COIL', 'REGCOIL_C2P_DIST_MIN', 'REGCOIL_LAMBDA')):
                print('Scan over ' + my_scan_item + ' with ' + str(my_scan_count) + ' steps')
                # Loop over number of variations
                for ii in range(0, my_scan_count):
                    # Make the 'edits'
                    self.OPTIMUM_TARGETS_PARAMS[my_scan_item]['TARGET'].set(self.OPTIMUM_SCAN_PARAMS[my_scan_item]['values'][ii])
                    self.OPTIMUM_TARGETS_PARAMS[my_scan_item]['SIGMA'].set(self.OPTIMUM_SCAN_PARAMS[my_scan_item]['sigmas'][ii])
                    # Modify the foldername_in
                    if (my_scan_item == 'REGCOIL_MAX_K'):
                        fext = '_RMK'
                    elif (my_scan_item == 'REGCOIL_RMS_K'):
                        fext = '_RRK'
                    elif (my_scan_item == 'REGCOIL_CHI2_B'):
                        fext = '_RC2B'
                    elif (my_scan_item == 'REGCOIL_BNORMAL_TOTAL'):
                        fext = '_RCBNT'
                    elif (my_scan_item == 'REGCOIL_MAX_BNORMAL'):
                        fext = '_RCMBN'
                    elif (my_scan_item == 'REGCOIL_VOLUME_COIL'):
                        fext = '_RVC'
                    elif (my_scan_item == 'REGCOIL_C2P_DIST_MIN'):
                        fext = '_RC2P'
                    elif (my_scan_item == 'REGCOIL_LAMBDA'):
                        fext = '_RL'
                    the_foldername_out = the_foldername_in + fext + str(ii)
                    self.recursive_scan_variation(variation_list_in, 
                                     variation_count_in, 
                                     the_foldername_out,
                                     the_filename_in)


            if (my_scan_item in('REGCOIL_FOURIER_SPECTRUM')):
                print('Scan over: ' +
                       str(self.OPTIMUM_SCAN_PARAMS['REGCOIL_FOURIER_SPECTRUM']['Selections']))
                # Loop over number of variations
                rfs_count = 0
                for rfs_selection in self.OPTIMUM_SCAN_PARAMS['REGCOIL_FOURIER_SPECTRUM']['Selections']:
                    # Make the 'edits'
                    if (self.OPTIMUM_SCAN_PARAMS['CLEAR_VARS_B4_REGCOIL']['Enabled'].get() is True):
                        self.clear_variables()

                    if (rfs_selection == 'REGCOIL_2x2'):
                        self.append_variables_regcoil_2x2()
                    elif (rfs_selection == 'REGCOIL_4x4'):
                        self.append_variables_regcoil_4x4()
                    elif (rfs_selection == 'REGCOIL_6x6'):
                        self.append_variables_regcoil_6x6()
                    elif (rfs_selection == 'REGCOIL_8x8'):
                        self.append_variables_regcoil_8x8()
                    elif (rfs_selection == 'REGCOIL_10x10'):
                        self.append_variables_regcoil_10x10()
                    elif (rfs_selection == 'REGCOIL_12x12'):
                        self.append_variables_regcoil_12x12()
                    elif (rfs_selection == 'REGCOIL_16x16'):
                        self.append_variables_regcoil_16x16()
                    elif (rfs_selection == 'REGCOIL_24x16'):
                        self.append_variables_regcoil_24x16()
                    else:
                        print('<----Something went wrong when parsing rfs_selection in recursive_scan_variation')
                    
                    # Modify the foldername_in
                    the_foldername_out = the_foldername_in + '_RFS' + str(rfs_count)
                    rfs_count += 1
                    
                    # append bounds, if desired
                    if (self.OPTIMUM_SCAN_PARAMS['AUTOGEN_REGCOIL_BOUNDS']['Enabled'].get() is True):
                        self.append_regcoil_bounds()

                    # append bounds, if desired
                    if (self.OPTIMUM_SCAN_PARAMS['AUTOGEN_REGCOIL_D']['Enabled'].get() is True):
                        self.append_regcoil_d()
                    
                    # Call resursively with reduced list
                    self.recursive_scan_variation(variation_list_in, 
                                     variation_count_in, 
                                     the_foldername_out,
                                     the_filename_in)

            if (my_scan_item in('FACTOR', 'EPSFCN')):
                self.OPTIMUM_SCAN_PARAMS[my_scan_item]['values']

                print('Scan over ' + my_scan_item + ' with ' + str(my_scan_count) + ' steps')
                # Loop over number of variations
                for ii in range(0, my_scan_count):
                    # Make the 'edits'
                    self.OPTIMUM_PARAMS[my_scan_item].set(self.OPTIMUM_SCAN_PARAMS[my_scan_item]['values'][ii])
                    # Modify the foldername_in
                    if (my_scan_item == 'FACTOR'):
                        fext = '_F'
                    elif (my_scan_item == 'EPSFCN'):
                        fext = '_E'
                    the_foldername_out = the_foldername_in + fext + str(ii)
                    self.recursive_scan_variation(variation_list_in, 
                                     variation_count_in, 
                                     the_foldername_out,
                                     the_filename_in)

            if (my_scan_item in('BOUNDARY_SPECTRUM')):
                print('Scan over: ' +
                       str(self.OPTIMUM_SCAN_PARAMS['BOUNDARY_SPECTRUM']['Selections']))
                # Loop over number of variations
                lcfs_count = 0
                for lcfs_selection in self.OPTIMUM_SCAN_PARAMS['BOUNDARY_SPECTRUM']['Selections']:
                    # Make the 'edits'
                    if (self.OPTIMUM_SCAN_PARAMS['CLEAR_VARS_B4_BOUNDARY']['Enabled'].get() is True):
                        self.clear_variables()
                    self.clear_variables()

                    if (lcfs_selection == 'BOUND_2x2'):
                        print('<----2')
                        self.append_variables_boundary_2x2()
                    elif (lcfs_selection == 'BOUND_4x4'):
                        print('<----4')
                        self.append_variables_boundary_4x4()
                    elif (lcfs_selection == 'BOUND_6x6'):
                        print('<----6')
                        self.append_variables_boundary_6x6()
                    elif (lcfs_selection == 'BOUND_8x8'):
                        print('<----8')
                        self.append_variables_boundary_8x8()
                    else:
                        print('<----Something went wrong when parsing lcfs_selection in recursive_scan_variation')
                    
                    # Modify the foldername_in
                    the_foldername_out = the_foldername_in + '_LCFS' + str(lcfs_count)
                    lcfs_count += 1
                    
                    # append bounds, if desired
                    if (self.OPTIMUM_SCAN_PARAMS['AUTOGEN_BOUNDARY_LIMITS']['Enabled'].get() is True):
                        self.append_boundary_bounds()

                    # append D's, if desired
                    if (self.OPTIMUM_SCAN_PARAMS['AUTOGEN_BOUNDARY_D']['Enabled'].get() is True):
                        self.append_boundary_d()

                    # Call resursively with reduced list
                    self.recursive_scan_variation(variation_list_in, 
                                     variation_count_in, 
                                     the_foldername_out,
                                     the_filename_in)


            variation_list_in.append(my_scan_item)
            variation_count_in.append(my_scan_count)
     
        
    def append_regcoil_bounds(self):
        # scan over input nescin file and assign "reasonable" bounds        
        source_file = self.FILESETC['NESCIN_FILEIN'].get()
        nescin_file = open(source_file)
        match_pattern = "------ Current Surface"
        
        # skip to the line that matches the pattern
        next_line = nescin_file.readline()
        while (next_line.find(match_pattern) < 0):
            next_line = nescin_file.readline()
        
        # read in # of Fourier modes
        next_line = nescin_file.readline()
        next_line = nescin_file.readline()
        num_fourier_modes = int(next_line)
        
        # read in two text lines that don't contain useful information
        next_line = nescin_file.readline()
        next_line = nescin_file.readline()
        
        # read in components. 
        m = []
        n = []
        rmnc_coil = []
        zmns_coil = []
        rmns_coil = []
        zmnc_coil = []
        
        for ii in range(0,num_fourier_modes):
            next_line = nescin_file.readline()
            (newm, newn, newrmnc, newzmns, newrmns, newzmnc) = next_line.split()
            m.append(int(newm))
            n.append(int(newn))
            rmnc_coil.append(float(newrmnc))
            zmns_coil.append(float(newzmns))
            rmns_coil.append(float(newrmns))
            zmnc_coil.append(float(newzmnc))

        # be nice and close the file            
        nescin_file.close()

        # Generate bounds.
        bounds_text = ''
        
        for ii in range(0,num_fourier_modes):
            #print('<----ii: ' + str(ii) + ', rmnc: ' + str(abs(rmnc_coil[ii])))
            if (abs(rmnc_coil[ii]) < 0.0005):
                r_bounds_min = '  REGCOIL_RCWS_RBOUND_C_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.005'
                r_bounds_max = '  REGCOIL_RCWS_RBOUND_C_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.005'
            elif (abs(rmnc_coil[ii]) < 0.005):
                r_bounds_min = '  REGCOIL_RCWS_RBOUND_C_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.02'
                r_bounds_max = '  REGCOIL_RCWS_RBOUND_C_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.02'
            elif (abs(rmnc_coil[ii]) < 0.02):
                r_bounds_min = '  REGCOIL_RCWS_RBOUND_C_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.1'
                r_bounds_max = '  REGCOIL_RCWS_RBOUND_C_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.1'
            elif (abs(rmnc_coil[ii]) < 0.05):
                r_bounds_min = '  REGCOIL_RCWS_RBOUND_C_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.2'
                r_bounds_max = '  REGCOIL_RCWS_RBOUND_C_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.2'
            elif (abs(rmnc_coil[ii]) < 0.30):
                r_bounds_min = '  REGCOIL_RCWS_RBOUND_C_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.6'
                r_bounds_max = '  REGCOIL_RCWS_RBOUND_C_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.6'
            elif (abs(rmnc_coil[ii]) < 0.65):
                r_bounds_min = '  REGCOIL_RCWS_RBOUND_C_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.8'
                r_bounds_max = '  REGCOIL_RCWS_RBOUND_C_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.8'
            else:
                #print('m = ' + str(m[ii]) + ', n = ' + str(n[ii]))
                if ((m[ii] == 0) and (n[ii] == 0)):
                    r_min00 = 0.75 * rmnc_coil[ii]
                    r_max00 = 1.25 * rmnc_coil[ii]
                    r_bounds_min = '  REGCOIL_RCWS_RBOUND_C_MIN(' + str(m[ii]) + \
                                   ', ' + str(n[ii]) + ') = ' + str(r_min00)
                    r_bounds_max = '  REGCOIL_RCWS_RBOUND_C_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = ' + str(r_max00)
                else:
                    print('<----Check the nescin file.  Large RMNC Mode? (m,n) = (' + str(m[ii]) + ', ' + str(n[ii]) + ')')

            if (abs(zmns_coil[ii]) < 0.0005):
                z_bounds_min = '  REGCOIL_RCWS_ZBOUND_S_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.005'
                z_bounds_max = '  REGCOIL_RCWS_ZBOUND_S_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.005'
            elif (abs(zmns_coil[ii]) < 0.005):
                z_bounds_min = '  REGCOIL_RCWS_ZBOUND_S_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.02'
                z_bounds_max = '  REGCOIL_RCWS_ZBOUND_S_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.02'
            elif (abs(zmns_coil[ii]) < 0.02):
                z_bounds_min = '  REGCOIL_RCWS_ZBOUND_S_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.1'
                z_bounds_max = '  REGCOIL_RCWS_ZBOUND_S_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.1'
            elif (abs(zmns_coil[ii]) < 0.05):
                z_bounds_min = '  REGCOIL_RCWS_ZBOUND_S_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.2'
                z_bounds_max = '  REGCOIL_RCWS_ZBOUND_S_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.2'
            elif (abs(zmns_coil[ii]) < 0.30):
                z_bounds_min = '  REGCOIL_RCWS_ZBOUND_S_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.6'
                z_bounds_max = '  REGCOIL_RCWS_ZBOUND_S_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.6'
            elif (abs(zmns_coil[ii]) < 0.65):
                z_bounds_min = '  REGCOIL_RCWS_ZBOUND_S_MIN(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = -0.8'
                z_bounds_max = '  REGCOIL_RCWS_ZBOUND_S_MAX(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.8'
            else:
                print('<----Check the nescin file.  Large ZMNS Mode? (m,n) = (' + str(m[ii]) + ', ' + str(n[ii]) + ')')
                            
            bounds_text += r_bounds_min + '\n' + r_bounds_max + '\n'
            if ( (m[ii] == 0) and (n[ii] == 0) ):
                # (0,0) mode
                pass
            else:
                bounds_text += z_bounds_min + '\n' + z_bounds_max + '\n'
                              
        # write the bounds to the end of the appropriate variable
        new_text = tk.StringVar()
        new_text.set(bounds_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

        
    def append_regcoil_d(self):
        # scan over input nescin file and assign "reasonable" bounds        
        source_file = self.FILESETC['NESCIN_FILEIN'].get()
        nescin_file = open(source_file)
        match_pattern = "------ Current Surface"
        
        # skip to the line that matches the pattern
        next_line = nescin_file.readline()
        while (next_line.find(match_pattern) < 0):
            next_line = nescin_file.readline()
        
        # read in # of Fourier modes
        next_line = nescin_file.readline()
        next_line = nescin_file.readline()
        num_fourier_modes = int(next_line)
        
        # read in two text lines that don't contain useful information
        next_line = nescin_file.readline()
        next_line = nescin_file.readline()
        
        # read in components. 
        m = []
        n = []
        rmnc_coil = []
        zmns_coil = []
        rmns_coil = []
        zmnc_coil = []
        
        for ii in range(0,num_fourier_modes):
            next_line = nescin_file.readline()
            (newm, newn, newrmnc, newzmns, newrmns, newzmnc) = next_line.split()
            m.append(int(newm))
            n.append(int(newn))
            rmnc_coil.append(float(newrmnc))
            zmns_coil.append(float(newzmns))
            rmns_coil.append(float(newrmns))
            zmnc_coil.append(float(newzmnc))

        # be nice and close the file            
        nescin_file.close()

        # Generate bounds.
        dbounds_text = ''
        
        for ii in range(0,num_fourier_modes):
            #print('<----ii: ' + str(ii) + ', rmnc: ' + str(abs(rmnc_coil[ii])))
            if (abs(rmnc_coil[ii]) < 0.0005):
                dr_bounds = '  DREGCOIL_RCWS_RBOUND_C_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.005'
            elif (abs(rmnc_coil[ii]) < 0.005):
                dr_bounds = '  DREGCOIL_RCWS_RBOUND_C_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.02'
            elif (abs(rmnc_coil[ii]) < 0.02):
                dr_bounds = '  DREGCOIL_RCWS_RBOUND_C_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.1'
            elif (abs(rmnc_coil[ii]) < 0.05):
                dr_bounds = '  DREGCOIL_RCWS_RBOUND_C_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.2'
            elif (abs(rmnc_coil[ii]) < 0.30):
                dr_bounds = '  DREGCOIL_RCWS_RBOUND_C_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.6'
            elif (abs(rmnc_coil[ii]) < 0.65):
                dr_bounds = '  DREGCOIL_RCWS_RBOUND_C_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.8'
            else:
                #print('m = ' + str(m[ii]) + ', n = ' + str(n[ii]))
                if ((m[ii] == 0) and (n[ii] == 0)):
                    dr_bounds = '  DREGCOIL_RCWS_RBOUND_C_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = ' + str(rmnc_coil[ii])
                else:
                    print('<----Check the nescin file.  Large RMNC Mode? (m,n) = (' + str(m[ii]) + ', ' + str(n[ii]) + ')')

            if (abs(zmns_coil[ii]) < 0.0005):
                dz_bounds = '  DREGCOIL_RCWS_ZBOUND_S_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.005'
            elif (abs(zmns_coil[ii]) < 0.005):
                dz_bounds = '  DREGCOIL_RCWS_ZBOUND_S_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.02'
            elif (abs(zmns_coil[ii]) < 0.02):
                dz_bounds = '  DREGCOIL_RCWS_ZBOUND_S_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.1'
            elif (abs(zmns_coil[ii]) < 0.05):
                dz_bounds = '  DREGCOIL_RCWS_ZBOUND_S_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.2'
            elif (abs(zmns_coil[ii]) < 0.30):
                dz_bounds = '  DREGCOIL_RCWS_ZBOUND_S_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.6'
            elif (abs(zmns_coil[ii]) < 0.65):
                dz_bounds = '  DREGCOIL_RCWS_ZBOUND_S_OPT(' + str(m[ii]) + \
                               ', ' + str(n[ii]) + ') = 0.8'
            else:
                print('<----Check the nescin file.  Large ZMNS Mode? (m,n) = (' + str(m[ii]) + ', ' + str(n[ii]) + ')')
                            
            dbounds_text += dr_bounds + '\n'
            if ( (m[ii] == 0) and (n[ii] == 0) ):
                # (0,0) mode
                pass
            else:
                dbounds_text += dz_bounds + '\n'
                              
        # write the bounds to the end of the appropriate variable
        new_text = tk.StringVar()
        new_text.set(dbounds_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())
        

    def append_boundary_bounds(self):
        # scan over input  file and assign "reasonable" bounds        
        source_file = self.FILESETC['BOUNDARY_INIT_IN'].get()
        
        new_vmec_nml = f90nml.read(source_file)
        rbc = new_vmec_nml['input']['rbc']
        zbs = new_vmec_nml['input']['zbs']
        
        bounds_text = ''
        for nn in range (0, 9):
            for mm in range(-8,9):
                this_rbc = rbc[nn][mm+8]
                this_zbs = zbs[nn][mm+8]
                #print('  mm=', str(mm), ' nn=', str(nn), ', rbc(mm,nn) = ', str(this_rbc))
                #print('  mm=', str(mm), ' nn=', str(nn), ', zbs(mm,nn) = ', str(this_zbs))
                if (abs(this_rbc) < 0.0005):
                    rbc_min = '  RBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.005'
                    rbc_max = '  RBC_MAX(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.005'
                elif (abs(this_rbc) < 0.005):
                    rbc_min = '  RBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.02'
                    rbc_max = '  RBC_MAX(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.02'
                elif (abs(this_rbc) < 0.02):
                    rbc_min = '  RBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.1'
                    rbc_max = '  RBC_MAX(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.1'
                elif (abs(this_rbc) < 0.05):
                    rbc_min = '  RBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.2'
                    rbc_max = '  RBC_MAX(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.2'
                elif (abs(this_rbc) < 0.30):
                    rbc_min = '  RBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.6'
                    rbc_max = '  RBC_MAX(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.6'
                elif (abs(this_rbc) < 0.65):
                    rbc_min = '  RBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.8'
                    rbc_max = '  RBC_MAX(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.8'
                else:
                    #print('m = ' + str(m[ii]) + ', n = ' + str(n[ii]))
                    if ((mm == 0) and (nn == 0)):
                        r_min00 = 0.75 * this_rbc
                        r_max00 = 1.25 * this_rbc
                        rbc_min = '  RBC_MIN(' + str(mm) + \
                                       ', ' + str(nn) + ') = ' + str(r_min00)
                        rbc_max = '  RBC_MAX(' + str(mm) + \
                                       ', ' + str(nn) + ') = ' + str(r_max00)
                    else:
                        print('<----Check the input file.  Large RMNC Mode? (m,n) = (' + str(m[ii]) + ', ' + str(nn) + ')')
 
                if (abs(this_zbs) < 0.0005):
                    zbs_min = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.005'
                    zbs_max = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.005'
                elif (abs(this_zbs) < 0.005):
                    zbs_min = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.02'
                    zbs_max = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.02'
                elif (abs(this_zbs) < 0.02):
                    zbs_min = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.1'
                    zbs_max = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.1'
                elif (abs(this_zbs) < 0.05):
                    zbs_min = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.2'
                    zbs_max = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.2'
                elif (abs(this_zbs) < 0.30):
                    zbs_min = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.6'
                    zbs_max = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.6'
                elif (abs(this_zbs) < 0.65):
                    zbs_min = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = -0.8'
                    zbs_max = '  ZBC_MIN(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.8'
                else:
                    print('<----Check the nescin file.  Large ZMNS Mode? (m,n) = (' + str(mm) + ', ' + str(nn) + ')')
                                 
                bounds_text += rbc_min + '\n' + rbc_max + '\n'
                if ( (mm == 0) and (nn == 0) ):
                    # (0,0) mode
                    pass
                else:
                    bounds_text += zbs_min + '\n' + zbs_max + '\n'
                                   
        # write the bounds to the end of the appropriate variable
        new_text = tk.StringVar()
        new_text.set(bounds_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())

        
    def append_boundary_d(self):
        # scan over input nescin file and assign "reasonable" bounds      
        source_file = self.FILESETC['BOUNDARY_INIT_IN'].get()
        
        new_vmec_nml = f90nml.read(source_file)
        rbc = new_vmec_nml['input']['rbc']
        zbs = new_vmec_nml['input']['zbs']
        
        # Generate bounds.
        dbound_text = ''
        
        for nn in range (0, 9):
            for mm in range(-8,9):
                this_rbc = rbc[nn][mm+8]
                this_zbs = zbs[nn][mm+8]
                this_max = max([abs(this_zbs), abs(this_rbc)])

                #print('  mm=', str(mm), ' nn=', str(nn), ', rbc(mm,nn) = ', str(this_rbc))
                #print('  mm=', str(mm), ' nn=', str(nn), ', zbs(mm,nn) = ', str(this_zbs))
                if (this_max < 0.0005):
                    dbound = '  DBOUND(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.005'
                elif (this_max < 0.005):
                    dbound = '  DBOUND(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.02'
                elif (this_max < 0.020):
                    dbound = '  DBOUND(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.1'
                elif (this_max < 0.05):
                    dbound = '  DBOUND(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.2'
                elif (this_max < 0.30):
                    dbound = '  DBOUND(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.6'
                elif (this_max < 0.65):
                    dbound = '  DBOUND(' + str(mm) + \
                                   ', ' + str(nn) + ') = 0.8'
                else:
                    if ((mm == 0) and (nn == 0)):
                        dbound = '  DBOUND(' + str(mm) + \
                                   ', ' + str(nn) + ') = ' + str(this_max)
                    else:
                        print('<----Check the nescin file.  Large RMNC Mode? (m,n) = (' + str(mm) + ', ' + str(nn) + ')')
                             
                dbound_text += dbound + '\n'
                              
        # write the bounds to the end of the appropriate variable
        new_text = tk.StringVar()
        new_text.set(dbound_text)
        self.OPTIMUM_VARIABLE_LINES.insert('end', new_text.get())
        
    
    
    
    def write_indata_nml(self):
            new_vmec_nml = f90nml.Namelist()
            new_vmec_nml['indata'] = {}
            for this_item in ('ntor', 'mpol', 'nzeta', 'ntheta', 'ncurr', 'nstep', 'nvacskip', 'delt', 'tcon0',
                              'lforbal', 'lasym', 'nfp', 'phiedge'):
                try:
                    new_vmec_nml['indata'][this_item] = self.VMEC_RUN_PARAMS[str.upper(this_item)].get()
                except:
                    print('<----Did not find indata/: ' + this_item + '  Skipping')
                
            try:
                new_vmec_nml['indata']['mgrid_file'] = self.VMEC_RUN_PARAMS['MGRID_FILE'].get().strip()
            except:
                    print('<----Did not find: mgrid_file Skipping')
                    
            for this_item in ('ftol_array', 'EXTCUR'):
                try:
                    the_list = self.VMEC_RUN_PARAMS[str.upper(this_item)].get().replace(',', ' ').split()
                    the_numbers = []
                    for my_item in the_list:
                        the_numbers.append(float(my_item))
                    new_vmec_nml['indata'][this_item] = the_numbers
                except:         
                    print('<----Did not find indata/: ' + this_item + ' Skipping')

            for this_item in ('ns_array', 'niter_array'):
                try:
                    the_list = self.VMEC_RUN_PARAMS[str.upper(this_item)].get().replace(',', ' ').split()
                    the_numbers = []
                    for my_item in the_list:
                        the_numbers.append(int(my_item))
                    new_vmec_nml['indata'][this_item] = the_numbers
                except:         
                    print('<----Did not find indata/: ' + this_item + ' Skipping')
                    
            new_vmec_nml.write(self.my_tempfile1, force=True)

            # append the initial axis, boundary and coils data
            # to the end of the file

            # read the file and output the text right back to the file, skipping the
            # 'trailing backslash'
            input_file = open(self.my_tempfile1, 'r')
            file_text = input_file.readlines()
            input_file.close()
            input_file = open(self.my_tempfile1, 'w')
            for this_line in file_text:
                if (this_line != '/\n'):
                    # in most cases, just print the line
                    input_file.write(this_line)
                else:
                    # in the case of the trailing backslash, print a newline
                    input_file.write('\n')

            # now do the append of the Initial Position info (axis and boundary)
            input_file.write(self.VMEC_InitPosition.get(1.0, "end"))

            input_file.write('/\n')

            # be nice. close the file handle
            input_file.close()
            return

    
    def write_optimum_nml(self):
        # check for 'coordinate' option
            new_optimum_nml = f90nml.Namelist()
            new_optimum_nml['optimum'] = {}
            new_optimum_nml['optimum']['NFUNC_MAX'] = self.OPTIMUM_PARAMS['NFUNC_MAX'].get()
            new_optimum_nml['optimum']['EQUIL_TYPE'] = self.OPTIMUM_PARAMS['EQUIL_TYPE'].get()
            new_optimum_nml['optimum']['OPT_TYPE'] = self.OPTIMUM_PARAMS['OPT_TYPE'].get()
            new_optimum_nml['optimum']['BOOTCALC_TYPE'] = self.OPTIMUM_PARAMS['BOOTCALC_TYPE'].get()
            new_optimum_nml['optimum']['FTOL'] = self.OPTIMUM_PARAMS['FTOL'].get()
            new_optimum_nml['optimum']['GTOL'] = self.OPTIMUM_PARAMS['GTOL'].get()
            new_optimum_nml['optimum']['FACTOR'] = self.OPTIMUM_PARAMS['FACTOR'].get()
            new_optimum_nml['optimum']['EPSFCN'] = self.OPTIMUM_PARAMS['EPSFCN'].get()
            new_optimum_nml['optimum']['MODE'] = self.OPTIMUM_PARAMS['MODE'].get()
            new_optimum_nml['optimum']['LKEEP_MINS'] = self.OPTIMUM_PARAMS['LKEEP_MINS'].get()

            new_optimum_nml.write(self.my_tempfile2, force=True)

            # read the file and output the text right back to the file, skipping the
            # 'trailing backslash'
            input_file = open(self.my_tempfile2, 'r')
            file_text = input_file.readlines()
            input_file.close()
            input_file = open(self.my_tempfile2, 'w')
            for this_line in file_text:
                if (this_line != '/\n'):
                    # in most cases, just print the line
                    input_file.write(this_line)
                else:
                    # in the case of the trailing backslash, print a newline
                    input_file.write('\n')
            input_file.write(self.OPTIMUM_EXTRA_LINES.get(1.0, "end"))
            
            if (self.OPTIMUM_SCAN_PARAMS['COORDINATE_REGCOIL_TARGET_VALUE']['Enabled'].get() is True):
                if (self.REGCOIL_PARAMS['target_option'].get() == 'rms_K'):
                    target_value = self.OPTIMUM_TARGETS_PARAMS['REGCOIL_RMS_K']['TARGET'].get()
                elif (self.REGCOIL_PARAMS['target_option'].get() == 'max_K'):
                    target_value = self.OPTIMUM_TARGETS_PARAMS['REGCOIL_MAX_K']['TARGET'].get()
                else:
                    print('<----Something is amiss with target_option. Check spelling: rms_K or  max_K')
                next_line = '  REGCOIL_TARGET_VALUE = ' + str(target_value) + '\n'
                input_file.write(next_line)
            
            # now write targets
            
            # REGCOIL
            for this_key in ('REGCOIL_MAX_K', 'REGCOIL_RMS_K', 'REGCOIL_CHI2_B', 'REGCOIL_BNORMAL_TOTAL', 'REGCOIL_MAX_BNORMAL',
                             'REGCOIL_VOLUME_COIL', 'REGCOIL_C2P_DIST_MIN', 'REGCOIL_LAMBDA'):
                count = self.OPTIMUM_TARGETS_PARAMS[this_key]['COUNT'].get()
                target = self.OPTIMUM_TARGETS_PARAMS[this_key]['TARGET'].get()
                sigma = self.OPTIMUM_TARGETS_PARAMS[this_key]['SIGMA'].get()
                if (count == 0):
                    pass
                else:
                    if (count == 1):
                        paren_part = '(1)'
                        mult_part = '1*'
                    else:
                        paren_part = '(1:' + str(count) + ')'
                        mult_part = str(count) + '*'
                    next_line = ('target_' + this_key + paren_part + ' = ' +
                                mult_part + str(target)  + '\n' +
                                'sigma_' + this_key + paren_part + ' = ' +
                                mult_part + str(sigma)  + '\n')
                    input_file.write(next_line)

            # COBRAVMEC
            this_key = 'BALLOON'
            count = self.OPTIMUM_TARGETS_PARAMS[this_key]['COUNT'].get()
            target = self.OPTIMUM_TARGETS_PARAMS[this_key]['TARGET'].get()
            sigma = self.OPTIMUM_TARGETS_PARAMS[this_key]['SIGMA'].get()
            if (count <= 1):
                pass
            else:
                next_line = ('  BALLOON_THETA = ' + 
                             self.OPTIMUM_TARGETS_PARAMS[this_key]['BALLOON_THETA'].get()  +
                              '\n' +
                              '  BALLOON_ZETA = ' + 
                             self.OPTIMUM_TARGETS_PARAMS[this_key]['BALLOON_ZETA'].get()  +
                              '\n')
                input_file.write(next_line)

                paren_part = '(1)'
                mult_part = '1*'
                next_line = ('TARGET_' + this_key + paren_part + ' = ' +
                            mult_part + str(target)  + '\n' +
                            'SIGMA_' + this_key + paren_part + ' = ' +
                            mult_part + str(sigma)  + '\n')
                input_file.write(next_line)

                paren_part = '(2:' + str(count) + ')'
                mult_part = str(count-1) + '*'
                next_line = ('TARGET_' + this_key + paren_part + ' = ' +
                            mult_part + str(target)  + '\n' +
                            'SIGMA_' + this_key + paren_part + ' = ' +
                            mult_part + str(sigma)  + '\n')
                input_file.write(next_line)

            # BOOZ_XFORM
            next_line = ('  MBOZ = ' + 
                         str(self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD']['MBOZ'].get() ) +
                          '\n' +
                          '  NBOZ = ' + 
                         str(self.OPTIMUM_TARGETS_PARAMS['BOOZER_COORD']['NBOZ'].get() ) +
                          '\n')
            input_file.write(next_line)

            # MAGWELL
            this_key = 'MAGWELL'
            count = self.OPTIMUM_TARGETS_PARAMS[this_key]['COUNT'].get()
            target = self.OPTIMUM_TARGETS_PARAMS[this_key]['TARGET'].get()
            sigma = self.OPTIMUM_TARGETS_PARAMS[this_key]['SIGMA'].get()
            if (count == 0):
                pass
            else:
                if (count == 1):
                    paren_part = '(1)'
                    mult_part = '1*'
                else:
                    paren_part = '(2:' + str(count) + ')'
                    mult_part = str(count-1) + '*'
                next_line = ('target_' + this_key + paren_part + ' = ' +
                            mult_part + str(target)  + '\n' +
                            'sigma_' + this_key + paren_part + ' = ' +
                            mult_part + str(sigma)  + '\n')
                input_file.write(next_line)

            # HELICITY
            input_file.write(self.Helicity_Text.get(1.0, "end"))

            # NEO
            input_file.write(self.NEO_Text.get(1.0, "end"))

            # GAMMA_C
            input_file.write(self.GAMMA_C_Text.get(1.0, "end"))

            # Aspect Ratio
            input_file.write(self.ASPECT_Text.get(1.0, "end"))

            # Profiles
            if (self.BOOTSTRAP_PARAMS['PROFILES']['Enabled'] is True):
               input_file.write(self.Profiles_Text.get(1.0, "end"))

            # Write the variable lines to the files
            input_file.write(self.OPTIMUM_VARIABLE_LINES.get(1.0, "end"))  
                     
            input_file.write('/\n')

            # be nice. close the file handle
            input_file.close()
            return
    
    def write_regcoil_nml(self):
        
        # check for 'coordinate' option
        new_regcoil_nml = f90nml.Namelist()
        new_regcoil_nml['regcoil_nml'] = {}
        new_regcoil_nml['regcoil_nml']['nlambda'] = self.REGCOIL_PARAMS['nlambda'].get()
        new_regcoil_nml['regcoil_nml']['ntheta_plasma'] = self.REGCOIL_PARAMS['ntheta_plasma'].get()
        new_regcoil_nml['regcoil_nml']['nzeta_plasma'] = self.REGCOIL_PARAMS['nzeta_plasma'].get()
        new_regcoil_nml['regcoil_nml']['ntheta_coil'] = self.REGCOIL_PARAMS['ntheta_coil'].get()
        new_regcoil_nml['regcoil_nml']['nzeta_coil'] = self.REGCOIL_PARAMS['nzeta_coil'].get()
        new_regcoil_nml['regcoil_nml']['mpol_potential'] = self.REGCOIL_PARAMS['mpol_potential'].get()
        new_regcoil_nml['regcoil_nml']['ntor_potential'] = self.REGCOIL_PARAMS['ntor_potential'].get()
        new_regcoil_nml['regcoil_nml']['general_option'] = self.REGCOIL_PARAMS['general_option'].get()
        new_regcoil_nml['regcoil_nml']['geometry_option_plasma'] = self.REGCOIL_PARAMS['geometry_option_plasma'].get()
        new_regcoil_nml['regcoil_nml']['geometry_option_coil'] = self.REGCOIL_PARAMS['geometry_option_coil'].get()
        new_regcoil_nml['regcoil_nml']['target_option'] = self.REGCOIL_PARAMS['target_option'].get()
        new_regcoil_nml['regcoil_nml']['nescin_filename'] = self.REGCOIL_PARAMS['nescin_filename'].get()
        new_regcoil_nml['regcoil_nml']['net_poloidal_current_Amperes'] = self.REGCOIL_PARAMS['net_poloidal_current_Amperes'].get()
        new_regcoil_nml['regcoil_nml']['net_toroidal_current_Amperes'] = self.REGCOIL_PARAMS['net_toroidal_current_Amperes'].get()
        new_regcoil_nml['regcoil_nml']['symmetry_option'] = self.REGCOIL_PARAMS['symmetry_option'].get()
        new_regcoil_nml['regcoil_nml']['load_bnorm'] = self.REGCOIL_PARAMS['load_bnorm'].get()
        new_regcoil_nml['regcoil_nml']['bnorm_filename'] = self.REGCOIL_PARAMS['bnorm_filename'].get()

        new_regcoil_nml.write(self.my_tempfile3, force=True)

        return

    def write_bootsj_nml(self):
        
        new_bootsj_nml = f90nml.Namelist()
        new_bootsj_nml['bootin'] = {}
        new_bootsj_nml['bootin']['mbuse'] = self.BOOTSTRAP_PARAMS['MBUSE'].get()
        new_bootsj_nml['bootin']['nbuse'] = self.BOOTSTRAP_PARAMS['NBUSE'].get()
        new_bootsj_nml['bootin']['zeff1'] = self.BOOTSTRAP_PARAMS['ZEFF1'].get()
        new_bootsj_nml['bootin']['dens0'] = self.BOOTSTRAP_PARAMS['DENS0'].get()
        new_bootsj_nml['bootin']['teti'] = self.BOOTSTRAP_PARAMS['TETI'].get()
        new_bootsj_nml['bootin']['tempres'] = self.BOOTSTRAP_PARAMS['TEMPRES'].get()
        new_bootsj_nml['bootin']['damp_bs'] = self.BOOTSTRAP_PARAMS['DAMP_BS'].get()
        new_bootsj_nml['bootin']['isymm0'] = self.BOOTSTRAP_PARAMS['ISYMM0'].get()

        for this_item in ('ate', 'ati'):
            try:
                the_list = self.BOOTSTRAP_PARAMS[str.upper(this_item)].get().replace(',', ' ').split()
                the_numbers = []
                for my_item in the_list:
                    the_numbers.append(float(my_item))
                new_bootsj_nml['bootin'][this_item] = the_numbers
            except:         
                print('<----Did not find: ' + this_item + ' Skipping')

        new_bootsj_nml.write(self.my_tempfile4, force=True)

        return

    def write_sfincs_nml(self):
        
        input_file = open(self.my_tempfile5, 'w')

        input_file.write(self.SFINCS_Text.get(1.0, "end"))

        input_file.close()
        return

    def create_files(self, folder_name_in, filename_in):
        output_directory = os.path.join('output', folder_name_in)
        
        if (os.path.isdir(output_directory)):
            print('<----Directory exists: ' + output_directory)
        else:
            print('<----Creating directory: ' + output_directory)
            os.mkdir(output_directory)
        
        # generate filename for the input.vmec file, and open it for writing
        input_fn_complete = os.path.join(output_directory, filename_in)
        input_file = open(input_fn_complete, 'w')
        
        # copy contents of temp files into 'input_file'.
        temp_file = open(self.my_tempfile1, 'r')
        
        file_text = temp_file.readlines()
        temp_file.close()
        for this_line in file_text:
            input_file.write(this_line)

        temp_file = open(self.my_tempfile2, 'r')
        
        file_text = temp_file.readlines()
        temp_file.close()
        for this_line in file_text:
            input_file.write(this_line)

        temp_file = open(self.my_tempfile3, 'r')
        
        file_text = temp_file.readlines()
        temp_file.close()
        for this_line in file_text:
            input_file.write(this_line)

        if (self.BOOTSTRAP_PARAMS['BOOTSJ']['Enabled'] is True):
           temp_file = open(self.my_tempfile4, 'r')
        
           file_text = temp_file.readlines()
           temp_file.close()
           for this_line in file_text:
               input_file.write(this_line)

        if (self.BOOTSTRAP_PARAMS['SFINCS']['Enabled'] is True):
           temp_file = open(self.my_tempfile5, 'r')
        
           file_text = temp_file.readlines()
           temp_file.close()
           for this_line in file_text:
               input_file.write(this_line)
               
        
        input_file.write('\n&END\n')
            
        # be nice, close the input file
        input_file.close()
            
        # delete temp files
        os.remove(self.my_tempfile1)
        os.remove(self.my_tempfile2)
        os.remove(self.my_tempfile3)
        if (self.BOOTSTRAP_PARAMS['BOOTSJ']['Enabled'] is True):
          os.remove(self.my_tempfile4)
        if (self.BOOTSTRAP_PARAMS['SFINCS']['Enabled'] is True):
          os.remove(self.my_tempfile5)

        # generate filename for the README file, and open it for writing
        readme_fn_complete = os.path.join(output_directory, 'README.txt')
        readme_file = open(readme_fn_complete, 'w')
        # Write out the 'README' text and close the file
        readme_file.write(self.README_LINES.get(1.0, "end"))
        readme_file.close()
             
        # generate filename for the CHTC submit file, and open it for writing
        chtc_fn_complete = os.path.join(output_directory, self.FILESETC['chtc'].get())
        chtc_file = open(chtc_fn_complete, 'w')
        # Write out the 'chtc_contents' text and close the file
        chtc_file.write(self.CHTC_CONTENTS.get(1.0, "end"))
        chtc_file.close()

        # generate filename for the stellopt_submit file, and open it for writing
        submit_fn_complete = os.path.join(output_directory, self.FILESETC['submit'].get())
        submit_file = open(submit_fn_complete, 'w')
        # Write out the 'stellopt_contents' text and close the file
        submit_file.write(self.SUBMIT_CONTENTS.get(1.0, "end"))
        submit_file.close()

        # copy additional required files
        if (self.FILESETC['NESCIN']['Enabled'].get() is True):
            source_file = self.FILESETC['NESCIN_FILEIN'].get()
            destination_file = os.path.join(output_directory, self.FILESETC['NESCIN_FILEOUT'].get())
            shutil.copyfile(source_file, destination_file)
                     
# This is where the magic happens
root = tk.Tk()
my_gui = stellgen(root)
root.mainloop()

