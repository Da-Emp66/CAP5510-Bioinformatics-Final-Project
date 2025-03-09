
import os
import shutil
from typing import Literal, Union
import pandas as pd
import requests
from tqdm import tqdm

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

class CASPTestSet:
    def __init__(
        self,
        filepath: str = "casp10_to_14_dataset.json",
        casp_version: Literal[
            "CASP10",
            "CASP11",
            "CASP12",
            "CASP13",
            "CASP14",
            "COMBINED",
        ] = "COMBINED",
        **download_kwargs,
    ):
        dataframes = []

        if not os.path.exists(filepath):

            # Download the dataset
            if casp_version == "COMBINED":
                for casp_version in ["CASP10", "CASP11", "CASP12", "CASP13", "CASP14"]:
                    dataframes.append(
                        self._download(
                            casp_version.lower(),
                            delete_temporary_folder=(False if casp_version != "CASP14" else True),
                            **download_kwargs
                        )
                    )
            else:
                dataframes.append(self._download(casp_version.lower(), **download_kwargs))
            
            # Concatenate all the datasets curated
            self.df = pd.concat(dataframes, ignore_index=True)

            self.filepath = filepath
            self.save(self.filepath)

        else:

            # Load from an existing file
            self.filepath = filepath
            self.df = self.load_from_file(self.filepath)

    def _download(
        self,
        target_folder: str,
        temporary_dir: str = "temporary_dir",
        delete_temporary_folder: bool = True
    ) -> pd.DataFrame:
        # DataFrame to write to
        write_df = pd.DataFrame(columns=[
            "target_id",
            "pdb_id",
            "real_fasta",
            "real_pdb",
        ])

        # Set the path to the CSV
        target_csv_path = f"CASP-Datasets/data/{target_folder}"

        # Set the temporary clone directory
        if len(temporary_dir) > 0:
            self.temporary_dir = temporary_dir
        else:
            self.temporary_dir = temporary_dir

        # Download the dataset
        if not os.path.exists(os.path.join(temporary_dir, target_csv_path)):
            os.makedirs(temporary_dir, exist_ok=False)
            os.system(f"cd {temporary_dir} && git clone https://github.com/Eryk96/CASP-Datasets.git && cd ..")

        # Read the target CSV
        read_df = pd.read_csv(os.path.join(temporary_dir, *target_csv_path.split("/"), "domain_summary.csv"))
        read_df_length = len(read_df.index)

        # Obtain the necessary information
        for _idx, row in tqdm(read_df.iterrows(), desc=target_folder, total=read_df_length):
            target_name = row["target"]
            pdb_name = row["pdb"]
            real_fasta = open(
                os.path.join(
                    temporary_dir,
                    *target_csv_path.split("/"),
                    "fasta",
                    f"{pdb_name}.fasta"
                )
            ).readlines()[1]
            real_pdb = download_real_pdb(protein_name=row["pdb"])
            new_row = pd.DataFrame({
                "target_id": [target_name],
                "pdb_id": [pdb_name],
                "real_fasta": [real_fasta],
                "real_pdb": [real_pdb],
            })
            write_df = pd.concat([write_df, new_row], ignore_index=True)

        # Delete the temporary files
        if delete_temporary_folder:
            shutil.rmtree(temporary_dir)

        return write_df

    def __len__(self):
        return len(self.df.index)

    def __getitem__(self, idx: int):
        return self.df.iloc[[idx]]
    
    def save(self, filepath: str):
        if filepath.endswith(".json"):
            self.df.to_json(filepath)
        elif filepath.endswith(".csv"):
            self.df.to_csv(filepath)
        elif filepath.endswith(".parquet"):
            self.df.to_parquet(filepath)
        else:
            raise NotImplementedError()
        
    def load_from_file(self, filepath: str):
        df = None

        if filepath.endswith(".json"):
            df = pd.read_json(filepath)
        elif filepath.endswith(".csv"):
            df = pd.read_csv(filepath)
        elif filepath.endswith(".parquet"):
            df = pd.read_parquet(filepath)
        else:
            raise NotImplementedError()
        
        return df
    


if __name__ == "__main__":
    dataset = CASPTestSet()
    print(dataset.df)
