import os
import pandas as pd
import sys
from qc.utils import FormatDataframeSamples
from qc.qc_pipeline import qc_pipeline
from tkinter import filedialog, Tk, Button, ttk, Entry, StringVar, messagebox, Label

def import_list_ms_runs():
    global df
    global filepath
    global userIonsFileLabel
    csv_file_path = filedialog.askopenfilename()
    print(csv_file_path)
    separator = "\t"
    if csv_file_path.endswith('.csv'):
        separator = ','
    df = pd.read_csv(csv_file_path, sep=separator)

    df = FormatDataframeSamples(df, basePathCsvRuns=os.path.dirname(csv_file_path))

    if len(df) == 0:
        return
    
    filepath.set(os.path.dirname(csv_file_path))
    if os.path.exists(filepath.get() + "/User-Ions.csv"):
        userIonsFileLabel.set("User ions file: DETECTED User-Ions.csv")
    else:
        userIonsFileLabel.set("User ions file: NONE")
    
    #Clear the treeview list items
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

def import_list_ms_runs_clipboard():
    global df
    global filepath
    global userIonsFileLabel
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
    #Clear the treeview list items
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
    # Take the path of the first run as the result path:
    filepath.set(df["MSRUNPATH"][0])
    if os.path.exists(filepath.get() + "/User-Ions.csv"):
        userIonsFileLabel.set("User ions file: DETECTED User-Ions.csv")
    else:
        userIonsFileLabel.set("User ions file: NONE")

def call_backend_process():
    global df
    global filepath
    if df is None or len(df) == 0:
        messagebox.showerror("Please import a list of MS runs.")
        return
    
    if os.path.exists(filepath.get() + "/config.toml"):
        qc_pipeline(df, filepath.get(), filepath.get() + "/config.toml")
    else:
        qc_pipeline(df, filepath.get())


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    root = Tk()
    root.geometry('800x600')
    root.title("IonToolPack v1 | PeakQC")
    root.iconbitmap(sys.executable)

    tree = ttk.Treeview(root, height=20, show="headings")
    filepath = StringVar()
    userIonsFileLabel = StringVar()
    userIonsFileLabel.set("User ions file: NONE")

    Button(root, text='Import list of MS runs (.csv, .txt)', command=import_list_ms_runs).pack()
    Button(root, text='Paste MS runs from clipboard', command=import_list_ms_runs_clipboard).pack()
    tree.pack(fill="both")
    
    Label(root, text="Output path:").pack() # add label for file path
    Entry(root, textvariable=filepath).pack() # add text box for file path
    Label(root, textvariable=userIonsFileLabel).pack() # add label for user ions file
    Button(root, text='Process', command=call_backend_process).pack()

    root.mainloop()
