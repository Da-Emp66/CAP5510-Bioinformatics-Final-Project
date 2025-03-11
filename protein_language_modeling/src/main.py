from src.data import CASPTestSet
from src.interfaces import ProteinPredictionTask
from src.models import ESM3Model
from src.utils import ProteinComparator


def main():
    model = ESM3Model()
    comparator = ProteinComparator()
    test_set = CASPTestSet()

    first_sample = test_set[0]
    
    resultant_pdb = model(
        ProteinPredictionTask.STRUCTURE_PREDICTION,
        protein=first_sample["real_fasta"],
    )

    resulting_alignment = comparator.compute_score_and_alignment(
        first_sample["real_pdb"],
        resultant_pdb,
    )[0]

    comparator.visualize_alignment(resulting_alignment)

if __name__ == "__main__":
    main()

