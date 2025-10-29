import os
import pandas as pd
import sys
from qc.utils import FormatDataframeSamples, export_mza_metadata
from mirador.image_time_vs_mz import generate_TimeVsMzImages
from mirador.image_time_vs_arrival_time import generate_TimeVsArrivalTimeImages
from tkinter import filedialog, Tk, Button, ttk, Entry, StringVar, Label, Scrollbar, Frame, Text, END
import numpy as np
import re
from gui_tabs import call_backend_peakqc, call_backend_tandemmatch, call_backend_peakquant, call_backend_mirador, call_backend_comparefeatures

class MainApplication:
    def __init__(self, root):
        # Shared GUI objects
        self.root = root
        self.treeMsRuns = None
        # Shared data/state, dfMsRuns and outputFolder are shared across several tabs
        self.dfMsRuns = None
        self.outputFolder = StringVar()

        # get default config settings:
        configFile = os.path.join(os.getcwd(),"default-config.toml")
        with open(configFile, "r", encoding="utf-8") as f:
            content = f.read()
        # Pattern: Match each section (e.g., ["Mirador"]) until next section (starting with [), or EOF
        pattern = r'(?s)^#---\s*Mirador\s*\n(.*?)(?=^#---\s|\Z)'
        match = re.search(pattern, content, re.MULTILINE)
        self.params_mirador = "Section start '#--- Mirador' not found."
        if match:
            self.params_mirador = match.group(1).strip() 

        pattern = r'(?s)^#---\s*PeakQC\s*\n(.*?)(?=^#---\s|\Z)'
        match = re.search(pattern, content, re.MULTILINE)
        self.params_peakqc = "Section start '#--- PeakQC' not found."
        if match:
            self.params_peakqc = match.group(1).strip()

        pattern = r'(?s)^#---\s*TandemMatch\s*\n(.*?)(?=^#---\s|\Z)'
        match = re.search(pattern, content, re.MULTILINE)
        self.params_tandemmatch = "Section start '#--- TandemMatch' not found."
        if match:
            self.params_tandemmatch = match.group(1).strip() 

        pattern = r'(?s)^#---\s*PeakQuant\s*\n(.*?)(?=^#---\s|\Z)'
        match = re.search(pattern, content, re.MULTILINE)
        self.params_peakquant = "Section start '#--- PeakQuant' not found."
        if match:
            self.params_peakquant = match.group(1).strip()

        pattern = r'(?s)^#---\s*CompareFeatures\s*\n(.*?)(?=^#---\s|\Z)'
        match = re.search(pattern, content, re.MULTILINE)
        self.params_comparefeatures = "Section start '#--- CompareFeatures' not found."
        if match:
            self.params_comparefeatures = match.group(1).strip()


    def import_list_ms_runs(self):
        csv_file_path = filedialog.askopenfilename()
        print(csv_file_path)
        separator = "\t"
        if csv_file_path.endswith('.csv'):
            separator = ','
        df = pd.read_csv(csv_file_path, sep=separator)
        df = FormatDataframeSamples(df, basePathCsvRuns=os.path.dirname(csv_file_path))
        if len(df) == 0:
            return
        self.dfMsRuns = df
        self.outputFolder.set(os.path.dirname(csv_file_path))
         
        # Clear and update the treeview list items
        tree = self.treeMsRuns
        for item in tree.get_children():
            tree.delete(item)
        tree["columns"] = list(df.columns)
        for col in df.columns:
            tree.heading(col, text=col)
        for idx, (_, row) in enumerate(df.iterrows()):
            if idx % 2 == 0:
                tree.insert("", "end", text=str(idx), values=list(row), tags=("evenrow",))
            else:
                tree.insert("", "end", text=str(idx), values=list(row))


    def import_list_ms_runs_clipboard(self):
        runs = []
        paths = []
        cb = root.clipboard_get()
        for item in cb.split('\n'):
            runx = os.path.basename(item)
            runs.append(runx)
            paths.append(item.removesuffix(runx))
        if len(runs) == 0:
            return
        
        df = pd.DataFrame({"MSRUN": runs, "MSRUNPATH": paths})
        df = FormatDataframeSamples(df)
        self.dfMsRuns = df
        # Take the path of the first run as the result path:
        self.outputFolder.set(os.path.dirname((df["MSRUNPATH"][0])))

        # Clear and update the treeview list items
        tree = self.treeMsRuns
        for item in tree.get_children():
            tree.delete(item)
        tree["columns"] = list(df.columns)
        for col in df.columns:
            tree.heading(col, text=col)
        for idx, (_, row) in enumerate(df.iterrows()):
            if idx % 2 == 0:
                tree.insert("", "end", text=str(idx), values=list(row), tags=("evenrow",))
            else:
                tree.insert("", "end", text=str(idx), values=list(row))


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    root = Tk()
    root.geometry('850x650')
    root.title("IonToolPack v1")
    root.iconbitmap(sys.executable)
    # Set background color of the entire window to white
    root.configure(bg='white')

    # Create the main application instance
    app = MainApplication(root)

    Button(root, text='Import list of MS runs (.csv, .txt)', command=app.import_list_ms_runs, bg='#f5f5f5').pack(pady=2)
    Button(root, text='Paste MS runs from clipboard', command=app.import_list_ms_runs_clipboard, bg='#f5f5f5').pack(pady=0)

    # Add table for MS runs, Treeview object with scrollbar (Frame):
    tree_frame = Frame(root, bg='white')
    tree_frame.pack(fill="both", expand=False)
    # Add vertical and horizontal scrollbars
    v_scrollbar = Scrollbar(tree_frame, orient="vertical")
    h_scrollbar = Scrollbar(tree_frame, orient="horizontal")
    v_scrollbar.pack(side="right", fill="y")
    h_scrollbar.pack(side="bottom", fill="x")
    tree = ttk.Treeview(tree_frame, height=6, show="headings", yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set, selectmode="extended") 
    tree.pack(fill="both", expand=True, padx=10)
    v_scrollbar.config(command=tree.yview)
    h_scrollbar.config(command=tree.xview)
    app.treeMsRuns = tree

    Label(root, text="Output folder:").pack() # add label for file path
    Entry(root, textvariable=app.outputFolder, width=50).pack(fill='x', padx=20) # add text box for file path

    # Create custom styles
    style = ttk.Style()
    style.theme_use('default')
    # Define normal and bold fonts
    style.configure("TNotebook.Tab", font=("Helvetica", 10))
    style.map("TNotebook.Tab",
          font=[("selected", ("Helvetica", 10, "bold"))])

    # ------------------------------------------------
    # Add tabs for tools, Notebook object:------------
    notebook = ttk.Notebook(root)
    notebook.pack(fill='both', expand=True, padx=20, pady=10)

    tab_mirador = ttk.Frame(notebook, style="Custom.TFrame")
    tab_peakqc = ttk.Frame(notebook, style="Custom.TFrame")
    tab_tandemmatch = ttk.Frame(notebook, style="Custom.TFrame")
    tab_peakquant = ttk.Frame(notebook, style="Custom.TFrame")
    tab_comparefeatures = ttk.Frame(notebook, style="Custom.TFrame")

    notebook.add(tab_mirador, text="Mirador")
    notebook.add(tab_peakqc, text="PeakQC")
    notebook.add(tab_tandemmatch, text="TandemMatch")
    notebook.add(tab_peakquant, text="PeakQuant")
    notebook.add(tab_comparefeatures, text="  |  CompareFeatures")
    
    # Mirador tab:
    # Add a label
    ttk.Label(tab_mirador, text="Parameters:").pack(pady=5)
    # Add the text box (multiline) with horizontal scrollbar and no softwrap
    tbox_params_mirador_frame = Frame(tab_mirador)
    tbox_params_mirador_frame.pack(padx=10, pady=5, fill="x")
    tbox_params_mirador = Text(tbox_params_mirador_frame, height=10, width=100, wrap='none')
    tbox_params_mirador.pack(side="left", fill="x", expand=True)
    tbox_params_mirador.insert("1.0", app.params_mirador)
    tbox_params_mirador_xscroll = Scrollbar(tbox_params_mirador_frame, orient="horizontal", command=tbox_params_mirador.xview)
    tbox_params_mirador_xscroll.pack(side="bottom", fill="x")
    tbox_params_mirador.config(xscrollcommand=tbox_params_mirador_xscroll.set)
    # Create a frame to hold the buttons
    button_frame_mirador = ttk.Frame(tab_mirador)
    button_frame_mirador.pack(pady=(5, 0))  
    
    # First row: All three display buttons in one row
    display_button_row_mirador = ttk.Frame(button_frame_mirador)
    display_button_row_mirador.pack(fill="x", pady=(0, 5))
    Button(display_button_row_mirador, text='Display MS1 & XIC', 
        command=lambda: call_backend_mirador(app.dfMsRuns, 
                                                 app.outputFolder.get(), 
                                                 tbox_params_mirador.get("1.0", END), 
                                                 mode="display_xics"), 
                                                 bg='#f5f5f5').pack(side="left", padx=(0, 10))
    Button(display_button_row_mirador, text='Display MS2', 
        command=lambda: call_backend_mirador(app.dfMsRuns, 
                                                 app.outputFolder.get(), 
                                                 tbox_params_mirador.get("1.0", END), 
                                                 mode="display_ms2"), 
                                                 bg='#f5f5f5').pack(side="left", padx=(0, 10))
    Button(display_button_row_mirador, text='Display XIM heatmap', 
        command=lambda: call_backend_mirador(app.dfMsRuns, 
                                                 app.outputFolder.get(), 
                                                 tbox_params_mirador.get("1.0", END), 
                                                 mode="display_xim_heatmap"), 
                                                 bg='#f5f5f5').pack(side="left")
    
    # Line separator
    separator_line_mirador = ttk.Frame(button_frame_mirador, height=40)
    separator_line_mirador.pack(fill="x", pady=1)
    
    # Second row: PDF Export buttons
    pdf_export_row_mirador = ttk.Frame(button_frame_mirador)
    pdf_export_row_mirador.pack(fill="x", pady=(0, 5))
    Button(pdf_export_row_mirador, text='Export MS1 & XIC (PDF)', 
        command=lambda: call_backend_mirador(app.dfMsRuns, 
                                                 app.outputFolder.get(), 
                                                 tbox_params_mirador.get("1.0", END),
                                                 mode="export_xics"), 
                                                 bg='#f5f5f5').pack(side="left", padx=(0, 10))
    Button(pdf_export_row_mirador, text='Export MS2 (PDF)', 
        command=lambda: call_backend_mirador(app.dfMsRuns, 
                                                 app.outputFolder.get(), 
                                                 tbox_params_mirador.get("1.0", END),
                                                 mode="export_ms2"), 
                                                 bg='#f5f5f5').pack(side="left", padx=(0, 10))
    Button(pdf_export_row_mirador, text='Export XIM heatmap (PDF)', 
        command=lambda: call_backend_mirador(app.dfMsRuns, 
                                                 app.outputFolder.get(), 
                                                 tbox_params_mirador.get("1.0", END),
                                                 mode="export_xim_heatmap"), 
                                                 bg='#f5f5f5').pack(side="left", padx=(0, 10))
    
    # Third row: Other Export buttons
    export_button_row_mirador = ttk.Frame(button_frame_mirador)
    export_button_row_mirador.pack(fill="x")
    Button(export_button_row_mirador, text='Export MZA metadata', 
        command=lambda: export_mza_metadata(app.outputFolder.get()), 
        bg='#f5f5f5').pack(side="left", padx=(0, 10))
    Button(export_button_row_mirador, text='Export time vs. m/z images', 
        command=lambda: generate_TimeVsMzImages(app.dfMsRuns, app.outputFolder.get()), 
        bg='#f5f5f5').pack(side="left", padx=(0, 10))
    Button(export_button_row_mirador, text='Export time vs. AT images', 
        command=lambda: generate_TimeVsArrivalTimeImages(app.dfMsRuns, app.outputFolder.get()), 
        bg='#f5f5f5').pack(side="left", padx=(0, 10))

    # PeakQC tab: 
    # Add a label
    ttk.Label(tab_peakqc, text="Parameters:").pack(pady=5)
    # Add the text box (multiline) for PeakQC parameters with horizontal scrollbar and no softwrap
    tbox_params_peakqc_frame = Frame(tab_peakqc)
    tbox_params_peakqc_frame.pack(padx=10, pady=5, fill="x")
    tbox_params_peakqc = Text(tbox_params_peakqc_frame, height=10, width=100, wrap='none')
    tbox_params_peakqc.pack(side="left", fill="x", expand=True)
    tbox_params_peakqc.insert("1.0", app.params_peakqc)
    tbox_params_peakqc_xscroll = Scrollbar(tbox_params_peakqc_frame, orient="horizontal", command=tbox_params_peakqc.xview)
    tbox_params_peakqc_xscroll.pack(side="bottom", fill="x")
    tbox_params_peakqc.config(xscrollcommand=tbox_params_peakqc_xscroll.set)
    # Create a frame to hold the buttons
    button_frame_peakqc = ttk.Frame(tab_peakqc)
    button_frame_peakqc.pack(side="bottom", pady=5)  
    # Add the 'Process' button with parameters from the textbox
    Button(button_frame_peakqc, text='Process', 
        command=lambda: call_backend_peakqc(app.dfMsRuns, 
                                                app.outputFolder.get(), 
                                                tbox_params_peakqc.get("1.0", END)), bg='#f5f5f5').pack(pady=2)


    # TandemMatch tab:
    # Add a label
    ttk.Label(tab_tandemmatch, text="Parameters:").pack(pady=5)
    # Add the text box (multiline) for TandemMatch parameters with horizontal scrollbar and no softwrap
    tbox_params_tandemmatch_frame = Frame(tab_tandemmatch)
    tbox_params_tandemmatch_frame.pack(padx=10, pady=5, fill="x")
    tbox_params_tandemmatch = Text(tbox_params_tandemmatch_frame, height=10, width=100, wrap='none')
    tbox_params_tandemmatch.pack(side="left", fill="x", expand=True)
    tbox_params_tandemmatch.insert("1.0", app.params_tandemmatch)
    tbox_params_tandemmatch_xscroll = Scrollbar(tbox_params_tandemmatch_frame, orient="horizontal", command=tbox_params_tandemmatch.xview)
    tbox_params_tandemmatch_xscroll.pack(side="bottom", fill="x")
    tbox_params_tandemmatch.config(xscrollcommand=tbox_params_tandemmatch_xscroll.set)
    # Create a frame to hold the buttons
    button_frame_tandemmatch = ttk.Frame(tab_tandemmatch)
    button_frame_tandemmatch.pack(side="bottom", pady=5)  
    # Add the 'Process' button
    Button(button_frame_tandemmatch, text='Process', 
        command=lambda: call_backend_tandemmatch(app.dfMsRuns, 
                                                    app.outputFolder.get(), 
                                                    tbox_params_tandemmatch.get("1.0", END)), 
                                                    bg='#f5f5f5').pack(pady=2)

    # PeakQuant tab:
    # Add a label
    ttk.Label(tab_peakquant, text="Parameters:").pack(pady=5)
    # Add the text box (multiline) for PeakQuant parameters with horizontal scrollbar and no softwrap
    tbox_params_peakquant_frame = Frame(tab_peakquant)
    tbox_params_peakquant_frame.pack(padx=10, pady=5, fill="x")
    tbox_params_peakquant = Text(tbox_params_peakquant_frame, height=10, width=100, wrap='none')
    tbox_params_peakquant.pack(side="left", fill="x", expand=True)
    tbox_params_peakquant.insert("1.0", app.params_peakquant)
    tbox_params_peakquant_xscroll = Scrollbar(tbox_params_peakquant_frame, orient="horizontal", command=tbox_params_peakquant.xview)
    tbox_params_peakquant_xscroll.pack(side="bottom", fill="x")
    tbox_params_peakquant.config(xscrollcommand=tbox_params_peakquant_xscroll.set)
    # Create a frame to hold the buttons
    button_frame_peakquant = ttk.Frame(tab_peakquant)
    button_frame_peakquant.pack(side="bottom", pady=5)  
    # Add the 'Process' button
    Button(button_frame_peakquant, text='Process', 
        command=lambda: call_backend_peakquant(app.dfMsRuns, 
                                                 app.outputFolder.get(), 
                                                 tbox_params_peakquant.get("1.0", END)), 
                                                 bg='#f5f5f5').pack(pady=2)
    
    # CompareFeatures tab:
    # Add a label
    ttk.Label(tab_comparefeatures, text="Parameters:").pack(pady=5)
    # Add the text box (multiline) for CompareFeatures parameters with horizontal scrollbar and no softwrap
    tbox_params_comparefeatures_frame = Frame(tab_comparefeatures)
    tbox_params_comparefeatures_frame.pack(padx=10, pady=5, fill="x")
    tbox_params_comparefeatures = Text(tbox_params_comparefeatures_frame, height=10, width=100, wrap='none')
    tbox_params_comparefeatures.pack(side="left", fill="x", expand=True)
    tbox_params_comparefeatures.insert("1.0", app.params_comparefeatures)
    tbox_params_comparefeatures_xscroll = Scrollbar(tbox_params_comparefeatures_frame, orient="horizontal", command=tbox_params_comparefeatures.xview)
    tbox_params_comparefeatures_xscroll.pack(side="bottom", fill="x")
    tbox_params_comparefeatures.config(xscrollcommand=tbox_params_comparefeatures_xscroll.set)
    # Create a frame to hold the buttons
    button_frame_comparefeatures = ttk.Frame(tab_comparefeatures)
    button_frame_comparefeatures.pack(side="bottom", pady=5)  
    # Add the 'Process' button
    Button(button_frame_comparefeatures, text='Process', 
        command=lambda: call_backend_comparefeatures(tbox_params_comparefeatures.get("1.0", END)), 
                                                       bg='#f5f5f5').pack(pady=2)

    root.mainloop()