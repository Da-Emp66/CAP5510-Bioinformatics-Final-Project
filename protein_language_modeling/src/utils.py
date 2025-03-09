from enum import Enum
import os
from typing import Optional, Union
import py3Dmol
import requests

class ProteinComparatorMethod(Enum):
    US_ALIGN = 0
    TM_ALIGN = 1
    ALL = 2

class ProteinComparator:
    def __init__(self, method: ProteinComparatorMethod):
        self.method = method

        if self.method == ProteinComparatorMethod.TM_ALIGN or self.method == ProteinComparatorMethod.ALL:
            self._ensure_script_exists(
                file_url=r"https://zhanggroup.org/TM-align/TMalign.cpp",
                compile_command="g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp",
            )
        
        if self.method == ProteinComparatorMethod.US_ALIGN or self.method == ProteinComparatorMethod.ALL:
            self._ensure_script_exists(
                file_url=r"https://zhanggroup.org/US-align/bin/module/USalign.cpp",
                compile_command="g++ -v -static -O3 -ffast-math -o USalign USalign.cpp",
            )

    def _ensure_script_exists(self, file_url: str, compile_command: Optional[str] = None):
        if not os.path.exists(file_url):
            self.download_file(file_url)
            if compile_command is not None:
                # Compile C++ file for algorithm
                os.system(compile_command)

    def download_file(self, url: str, file_to_write_to: Optional[str] = None):
        if file_to_write_to is None:
            filename = url.split("/")[-1]
        else:
            filename = file_to_write_to

        with open(filename, "w") as script_file:
            script_file.write(requests.get(url).text)
            script_file.close()

    def compute_alignment(self, pdb1: str, pdb2: str):
        pass

    def compute_score(self) -> float:
        pass


class MoleculeStructureVisualization:
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

