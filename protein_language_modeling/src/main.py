from data import CASPTestSet
from interfaces import ProteinPredictionTask
from models import ESM3Model
from utils import ProteinComparator


def main():
    # Instantiate components
    model = ESM3Model()
    comparator = ProteinComparator()
    test_set = CASPTestSet()

    # Grab the first data point from the test set
    first_sample = test_set[0]

    # run inference on the original protein FASTA
    print(first_sample["real_fasta"][0])
    resultant_pdb = model(
        ProteinPredictionTask.STRUCTURE_PREDICTION,
        protein=first_sample["real_fasta"][0][:-4],
        generation_config_kwargs={ "num_steps": 1, "temperature": 0.0 }
    )

    # Compute the optimal alignment and USalign/TMalign scores
    resulting_alignment = comparator.compute_score_and_alignment(
        resultant_pdb,
        first_sample["real_pdb"][0],
    )[0]

    # Print and visualize the protein alignment and scores
    print(f"TM-Score 1: {resulting_alignment.score1}")
    print(f"TM-Score 2: {resulting_alignment.score2}")
    comparator.visualize_alignment(resulting_alignment)

    # Compute the optimal alignment and USalign/TMalign scores
    resulting_alignment = comparator.compute_score_and_alignment(
        first_sample["real_pdb"][0],
        resultant_pdb,
    )[0]

    # Print and visualize the protein alignment and scores
    print(f"TM-Score 1: {resulting_alignment.score1}")
    print(f"TM-Score 2: {resulting_alignment.score2}")
    comparator.visualize_alignment(resulting_alignment)


if __name__ == "__main__":
    main()

