from io import StringIO
import os
from typing import Any, Dict, List, Optional, Set, Union
from dotenv import load_dotenv
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
import torch
import uuid
from interfaces import BaseProteinLanguageModel, ProteinPredictionReturnType, ProteinPredictionTask

load_dotenv()

class ESM3Model(BaseProteinLanguageModel):
    def __init__(
        self,
        model_id: str = "esm3-open",
        device: str = "cuda" if torch.cuda.is_available() else "cpu"
    ):
        login(token=os.getenv("HF_TOKEN"))
        self.model_id = model_id
        self.model: ESM3InferenceClient = ESM3.from_pretrained(model_id).to(device)
        
    def __call__(
        self,
        task: ProteinPredictionTask,
        protein: Union[str, Any],
        return_type: ProteinPredictionReturnType = ProteinPredictionReturnType.STRING,
        # Model-specific kwargs
        generation_config_kwargs: Optional[Dict[str, Any]] = { "num_steps": 8 },
    ) -> Union[str, Any]:
        if isinstance(protein, str):
            if task == ProteinPredictionTask.MASKED_SEQUENCE_COMPLETION or \
                task == ProteinPredictionTask.STRUCTURE_PREDICTION:
                protein = ESMProtein(sequence=protein)
            elif task == ProteinPredictionTask.INVERSE_FOLDING:
                temporary_file = StringIO(protein)
                protein = ESMProtein.from_pdb(temporary_file)
            elif task == ProteinPredictionTask.UNKNOWN:
                raise NotImplementedError()
        
        if task == ProteinPredictionTask.MASKED_SEQUENCE_COMPLETION:
            generation_config = GenerationConfig("sequence", temperature=0.7, **generation_config_kwargs)
            output: ESMProtein = self.model.generate(protein, generation_config)
            if return_type == ProteinPredictionReturnType.STRING:
                output: str = output.sequence
        elif task == ProteinPredictionTask.STRUCTURE_PREDICTION:
            protein.coordinates = None
            generation_config = GenerationConfig("structure", **generation_config_kwargs)
            output: ESMProtein = self.model.generate(protein, generation_config)
            if return_type == ProteinPredictionReturnType.STRING:
                output: str = output.to_pdb_string()
        elif task == ProteinPredictionTask.INVERSE_FOLDING:
            protein.sequence = None
            generation_config = GenerationConfig("sequence", **generation_config_kwargs)
            output: ESMProtein = self.model.generate(protein, generation_config)
            if return_type == ProteinPredictionReturnType.STRING:
                output: str = output.sequence

        return output
    
    def supported_tasks(self) -> Set[ProteinPredictionTask]:
        return set([
            ProteinPredictionTask.MASKED_SEQUENCE_COMPLETION,
            ProteinPredictionTask.STRUCTURE_PREDICTION,
            ProteinPredictionTask.INVERSE_FOLDING,
        ])

# class AlphaFold3Model(BaseProteinLanguageModel):
#     def __init__(self, use_server: bool = True):
#         self.use_server = use_server

#         if self.use_server:
#             pass
#         else:
#             raise NotImplementedError()

#     def __call__(
#         self,
#         task: ProteinPredictionTask,
#         protein: Union[str, Any],
#         return_type: ProteinPredictionReturnType = ProteinPredictionReturnType.STRING,
#         # Model-specific kwargs
#         model_seeds: List[int] = [1, 2],
#         protein_modifications: List[Any] = [],
#         unpaired_msa: Optional[str] = None,
#         paired_msa: Optional[str] = None,
#         protein_templates: List[Any] = [],
#         additional_job_kwargs: Dict[str, Any] = {},
#     ) -> Union[str, Any]:
        
#         if self.use_server:

#             job_name_id = str(uuid.uuid4())
#             protein_name_id = f"{job_name_id}-{str(uuid.uuid4())}"

#             model_input = {
#                 "name": job_name_id,
#                 "modelSeeds": model_seeds,
#                 "sequences": [
#                     {
#                         "protein": {
#                             "id": protein_name_id,
#                             "sequence": str(protein),
#                             "modifications": protein_modifications,
#                             "unpairedMsa": unpaired_msa,
#                             "pairedMsa": paired_msa,
#                             "templates": protein_templates,
#                         }
#                     },
#                 ],
#                 "dialect": "alphafold3",
#                 "version": 3,
#             }

#         else:

#             raise NotImplementedError()

#     def supported_tasks(self):
#         return set([])
    
if __name__ == "__main__":
    protein = "___________________________________________________DQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPP___________________________________________________________"
    model = ESM3Model() # "esm3_sm_open_v1"
    print("Step 1: Masked Sequence Completion")
    output = model(ProteinPredictionTask.MASKED_SEQUENCE_COMPLETION, protein)
    print(output)
    print("Step 2: Structure Prediction")
    output = model(ProteinPredictionTask.STRUCTURE_PREDICTION, output)
    print(output)
    print("Step 3: Inverse Folding Prediction")
    output = model(ProteinPredictionTask.INVERSE_FOLDING, output)
    print(output)

