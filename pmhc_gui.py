#!/usr/bin/env python3
"""Tkinter GUI for peptide-MHC-I binding prediction (tensorflow-free).

Uses the exported pmhc_models.npz weights and pure-numpy inference.
Run directly (python pmhc_gui.py) or as a bundled .exe.
"""
import os
import sys
import tkinter as tk
from tkinter import ttk

from pmhcbenchmark.numpy_model import load_numpy_models


def resource_path(name):
    # PyInstaller onefile extracts bundled data to sys._MEIPASS
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, name)


MODELS = load_numpy_models(resource_path("pmhc_models.npz"))
VALID_AA = "ARNDCQEGHILKMFPSTWYV"


class App:
    def __init__(self, root):
        root.title("pMHC-I Binding Predictor")
        root.geometry("460x240")
        root.resizable(False, False)

        ttk.Label(root, text="Peptide (9-mer):").grid(row=0, column=0, padx=12, pady=12, sticky="w")
        self.pep = ttk.Entry(root, width=22, font=("Consolas", 11))
        self.pep.grid(row=0, column=1, padx=12, pady=12)

        ttk.Label(root, text="Allele:").grid(row=1, column=0, padx=12, pady=6, sticky="w")
        self.allele = tk.StringVar(value="HLA-A*02:01")
        ttk.OptionMenu(root, self.allele, "HLA-A*02:01", *MODELS.keys()).grid(
            row=1, column=1, padx=12, pady=6, sticky="w")

        ttk.Button(root, text="Predict", command=self.predict).grid(
            row=2, column=0, columnspan=2, pady=14)

        self.result = tk.StringVar(value="")
        ttk.Label(root, textvariable=self.result, font=("Consolas", 11),
                  justify="left", anchor="w").grid(
            row=3, column=0, columnspan=2, padx=14, sticky="w")

        self.pep.bind("<Return>", lambda e: self.predict())

    def predict(self):
        pep = self.pep.get().strip().upper()
        if len(pep) != 9 or any(a not in VALID_AA for a in pep):
            self.result.set("Enter a valid 9-mer (20 standard amino acids).")
            return
        allele = self.allele.get()
        s = float(MODELS[allele].score([pep])[0])
        cls = "BINDER" if s >= 0.5 else "NON-BINDER"
        self.result.set(f"{pep}   {allele}\nscore = {s:.4f}   →   {cls}")


if __name__ == "__main__":
    root = tk.Tk()
    App(root)
    root.mainloop()
