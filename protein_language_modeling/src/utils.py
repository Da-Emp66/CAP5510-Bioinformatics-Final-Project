from dataclasses import dataclass
from enum import Enum
import os
import re
from typing import Any, Dict, List, Literal, Optional
import py3Dmol
import subprocess

import requests
from tempfile import NamedTemporaryFile
from esm.utils.structure.protein_chain import ProteinChain
# from tmscoring import TMscoring

class ProteinComparatorMethod(Enum):
    US_ALIGN = 0
    TM_ALIGN = 1
    # Root-mean-square Deviation Test
    RMSD = 2
    # Local Distance Difference Test
    LDDT = 3
    ALL = 4

class MoleculeStructureVisualization:
    def __init__(self, **view_kwargs):
        self._initialize(view_kwargs)

    def _initialize(self, view_kwargs: Dict[str, Any] = { "width": 500, "height": 500 }):
        self.initial_view_kwargs = view_kwargs
        self.view = py3Dmol.view(**view_kwargs)
        self.molecule_count = 0

    def add_molecule(
        self,
        pdb: str,
        color: Literal[
            "red",
            "blue",
            "yellow",
            "green",
            "spectrum"
        ] = "spectrum",
    ):
        self.view.addModel(pdb, "pdb")
        self.view.setStyle({"model": self.molecule_count}, {"cartoon": {"color": color}})
        self.molecule_count += 1

    def display(self):
        self.view.zoomTo()
        self.view.show()

    def reset(self):
        self._initialize(self.initial_view_kwargs)

@dataclass
class ProteinAlignment:
    method: ProteinComparatorMethod
    pdb1: str
    pdb2: str
    superimposed_pdb: str
    score1: float
    score2: float
    final_score: Optional[float]
    auxiliary: Optional[Any]

class ProteinComparator:
    def __init__(
        self,
        method: ProteinComparatorMethod = ProteinComparatorMethod.US_ALIGN,
        implementation: Literal["cpp", "py"] = "cpp",
        visualizer: MoleculeStructureVisualization = MoleculeStructureVisualization(),
    ):
        self.method = method
        self.implementation = implementation
        self.visualizer = visualizer

        if self.method == ProteinComparatorMethod.TM_ALIGN or self.method == ProteinComparatorMethod.ALL:
            self._ensure_script_exists(
                filename="TMalign",
                file_url=r"https://zhanggroup.org/TM-align/TMalign.cpp",
                compile_command="g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp",
            )
        
        if self.method == ProteinComparatorMethod.US_ALIGN or self.method == ProteinComparatorMethod.ALL:
            self._ensure_script_exists(
                filename="USalign",
                file_url=r"https://zhanggroup.org/US-align/bin/module/USalign.cpp",
                compile_command="g++ -v -static -O3 -ffast-math -o USalign USalign.cpp",
            )

    def _ensure_script_exists(self, filename: str, file_url: str, compile_command: Optional[str] = None):
        # Download the script at the given URL if it does not exists
        if not os.path.exists(filename):
            self.download_file(file_url)
            if compile_command is not None:
                # Compile C++ file for algorithm
                os.system(compile_command)

    def download_file(self, url: str, file_to_write_to: Optional[str] = None):
        # Assume the filename is the leaf of the protocol path
        # if not otherwise specified
        if file_to_write_to is None:
            filename = url.split("/")[-1]
        else:
            filename = file_to_write_to

        # Write the file to the current file system
        with open(filename, "w") as script_file:
            script_file.write(requests.get(url).text)
            script_file.close()

    def _run_cpp_executable(
        self,
        *file_args: str,
        run_command_template: str,
        output_filename: Optional[str] = None,
        suffix_to_add_if_output_file_not_found: str = ".pdb",
    ):
        # Place the contents of the args into temporary files
        # recognizable by the file system
        temporary_files = [NamedTemporaryFile("w+") for _ in file_args]
        for contents, file in zip(file_args, temporary_files):
            file.write(contents)

        file_arguments_for_command = ([file.name for file in temporary_files] + ([output_filename] if output_filename is not None else []))

        # Create the command
        command = (run_command_template.format(*file_arguments_for_command)).split(" ")

        # Run the script with the associated named files
        result = subprocess.run(
            command,
            capture_output=True,
            text=True
        )

        # Delete the temporary files
        for file in temporary_files:
            file.close()

        # If we cannot find the output file,
        # assume it was passed as a file without the extension
        if not os.path.exists(output_filename) and "." not in output_filename:
            output_filename += "." + suffix_to_add_if_output_file_not_found.strip(".")
        
        output_file_contents = None

        # Obtain the content of the output file, if it exists
        if output_filename is not None:
            with open(output_filename, "r") as output_file:
                output_file_contents = output_file.read()
                output_file.close()

            # Remove the output_filename now that we know
            # the contents in-memory
            os.remove(output_filename)

        # Return the superimposed PDB or other output file
        # and the output to stdout
        return (output_file_contents, result.stdout)

    def _run_score_alignment_algorithms(self, pdb1: str, pdb2: str) -> List[ProteinAlignment]:
        results = []

        # Compute TM-Alignment and Scores
        if self.method == ProteinComparatorMethod.TM_ALIGN or self.method == ProteinComparatorMethod.ALL:
            superimposed_pdb_string, printed_stdout = self._run_cpp_executable(
                pdb1,
                pdb2,
                run_command_template="./TMalign -mm 1 -ter 0 {} {} -o {}",
                output_filename="superimposed",
            )

            tm_scores = re.findall(r"TM-score\=\s*((?:\d|\.)+)", printed_stdout)

            results.append(
                ProteinAlignment(
                    method=ProteinComparatorMethod.TM_ALIGN,
                    pdb1=pdb1,
                    pdb2=pdb2,
                    superimposed_pdb=superimposed_pdb_string,
                    score1=tm_scores[0],
                    score2=tm_scores[1],
                    final_score=tm_scores[0],
                    auxiliary=printed_stdout,
                )
            )

        # Compute US-Alignment and Scores
        if self.method == ProteinComparatorMethod.US_ALIGN or self.method == ProteinComparatorMethod.ALL:
            superimposed_pdb_string, printed_stdout = self._run_cpp_executable(
                pdb1,
                pdb2,
                run_command_template="./USalign -mm 1 -ter 0 {} {} -o {}",
                output_filename="superimposed",
            )

            tm_scores = re.findall(r"TM-score\=\s*((?:\d|\.)+)", printed_stdout)

            results.append(
                ProteinAlignment(
                    method=ProteinComparatorMethod.US_ALIGN,
                    pdb1=pdb1,
                    pdb2=pdb2,
                    superimposed_pdb=superimposed_pdb_string,
                    score1=tm_scores[0],
                    score2=tm_scores[1],
                    final_score=tm_scores[0],
                    auxiliary=printed_stdout,
                )
            )

        # Compute Root-mean-square Deviation Scores
        if self.method == ProteinComparatorMethod.RMSD or self.method == ProteinComparatorMethod.ALL:
            pdb_file_1, pdb_file_2 =  (NamedTemporaryFile("w+"), NamedTemporaryFile("w+"))
            pdb_file_1.write(pdb1)
            pdb_file_2.write(pdb2)
            protein_chain_1, protein_chain_2 = (ProteinChain.from_pdb(pdb_file_1), ProteinChain.from_pdb(pdb_file_2))
            
            rmsd_1 = protein_chain_1.rmsd(protein_chain_2)
            rmsd_2 = protein_chain_2.rmsd(protein_chain_1)

            results.append(
                ProteinAlignment(
                    method=ProteinComparatorMethod.RMSD,
                    pdb1=pdb1,
                    pdb2=pdb2,
                    superimposed_pdb=None,
                    score1=rmsd_1,
                    score2=rmsd_2,
                    final_score=rmsd_1,
                    auxiliary=None,
                )
            )

        # Compute Local Distance Difference Test Scores
        if self.method == ProteinComparatorMethod.LDDT or self.method == ProteinComparatorMethod.ALL:
            pdb_file_1, pdb_file_2 =  (NamedTemporaryFile("w+"), NamedTemporaryFile("w+"))
            pdb_file_1.write(pdb1)
            pdb_file_2.write(pdb2)
            protein_chain_1, protein_chain_2 = (ProteinChain.from_pdb(pdb_file_1), ProteinChain.from_pdb(pdb_file_2))

            lddt_1 = protein_chain_1.lddt_ca(protein_chain_2, per_residue=False)
            lddt_2 = protein_chain_2.lddt_ca(protein_chain_1, per_residue=False)

            results.append(
                ProteinAlignment(
                    method=ProteinComparatorMethod.LDDT,
                    pdb1=pdb1,
                    pdb2=pdb2,
                    superimposed_pdb=None,
                    score1=lddt_1,
                    score2=lddt_2,
                    final_score=lddt_1,
                    auxiliary=None,
                )
            )

        return results

    def compute_score_and_alignment(self, pdb1: str, pdb2: str) -> List[ProteinAlignment]:
        """NOTE: The first PDB, `pdb1`, should always be the predicted PDB string,
        and `pdb2` should be the ground truth PDB string.

        Args:
            pdb1 (str): Predicted protein PDB representation.
            pdb2 (str): Ground truth protein PDB representation.

        Returns:
            List[ProteinAlignment]: All protein alignments and scoring metrics.
        """
        results = self._run_score_alignment_algorithms(pdb1, pdb2)
        return results
    
    def visualize_alignment(self, protein_alignment: ProteinAlignment, reference: Literal["pdb1", "pdb2"] = "pdb2", color1: str = "red", color2: str = "blue"):
        if reference == "pdb1":
            self.visualizer.add_molecule(protein_alignment.pdb1, color=color1)
        elif reference == "pdb2":
            self.visualizer.add_molecule(protein_alignment.pdb2, color=color1)

        self.visualizer.add_molecule(protein_alignment.superimposed_pdb, color=color2)
        self.visualizer.display()
        self.visualizer.reset()

    def visualize_comparison(self, pdb1: str, pdb2: str, color1: str = "red", color2: str = "blue"):
        self.visualizer.add_molecule(pdb1, color=color1)
        self.visualizer.add_molecule(pdb2, color=color2)
        self.visualizer.display()
        self.visualizer.reset()
