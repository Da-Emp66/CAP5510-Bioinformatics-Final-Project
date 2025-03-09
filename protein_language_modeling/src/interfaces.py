import abc
from enum import Enum
from typing import Any, Set, Union

class ProteinPredictionTask(Enum):
    MASKED_SEQUENCE_COMPLETION = 0
    """Predicting the masked elements of a protein sequence."""
    STRUCTURE_PREDICTION = 1
    """Predicting the three-dimensional coordinates of a protein sequence, usually in PDB format."""
    INVERSE_FOLDING = 2
    """Predicting the sequence of a protein based on its three-dimensional atomic structure."""
    UNKNOWN = 3
    """A prediction task determined by external arguments"""

class ProteinPredictionReturnType(Enum):
    STRING = str
    DEFAULT = Any
    UNKNOWN = Any

class BaseProteinLanguageModel(metaclass=abc.ABCMeta):
    def __init__(self):
        raise NotImplementedError
    
    def __call__(
        self,
        task: ProteinPredictionTask,
        protein: Union[str, Any],
        return_type: ProteinPredictionReturnType = ProteinPredictionReturnType.STRING,
        *args,
        **kwargs,
    ) -> Union[str, Any]:
        raise NotImplementedError
    
    def supported_tasks(self) -> Set[ProteinPredictionTask]:
        raise NotImplementedError
    