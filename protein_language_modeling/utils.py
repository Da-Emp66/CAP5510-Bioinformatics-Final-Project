import os
from typing import Union
import py3Dmol
import requests


class Visualization:
    def __init__(self):
        self.view = py3Dmol.view(width=500, height=500)

    def construct(self, pdb: str):
        self.view.addModel(pdb, "pdb")
        # Set the style of the protein chain
        self.view.setStyle({"cartoon": {"color": "spectrum"}})
        # Zoom in on the protein chain
        self.view.zoomTo()

    def display(self):
        # Display the protein chain
        self.view.show()


def download_real_pdb(protein_name: str) -> Union[str, None]:
    # Set the PDB to null for if it does not exist
    pdb = None

    # Download the PDB from RCSB (U.S. data center for Protein Data Bank [PDB])
    response = requests.get(f"https://files.rcsb.org/download/{protein_name}.pdb")

    # Set the PDB to the response PDB if the response code <= 399
    if response.ok:
        pdb = response.text

    # Return the PDB string
    return pdb

