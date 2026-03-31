import os
import sys
import pickle
from inspect import signature

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from rdkit import Chem, RDLogger
from moses.metrics.metrics import get_all_metrics

import datasets
from configs.vp_qm9_cdgs import get_config


ALLOWED_ATOMS = {"C", "N", "O", "F"}
MAX_ATOMS = 9

MAX_VALENCE = {
    "C": 4,
    "N": 4,   
    "O": 2,
    "F": 1,
}


def passes_atom_set(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ALLOWED_ATOMS:
            return False
    return True


def is_single_fragment(mol):
    return len(Chem.GetMolFrags(mol)) == 1


def passes_size_limit(mol, max_atoms=MAX_ATOMS):
    return mol.GetNumAtoms() <= max_atoms


def passes_valence_rules(mol):
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in MAX_VALENCE:
            return False

        explicit_valence = atom.GetExplicitValence()
        total_valence = atom.GetTotalValence()

        if explicit_valence > MAX_VALENCE[symbol]:
            return False
        if total_valence > MAX_VALENCE[symbol]:
            return False

    return True


def sanitize_and_filter(smiles_list):
    kept = []

    stats = {
        "input_count": len(smiles_list),
        "invalid_smiles": 0,
        "atom_set_fail": 0,
        "fragment_fail": 0,
        "size_fail": 0,
        "valence_fail": 0,
        "kept_count": 0,
    }

    for s in smiles_list:
        if s is None:
            stats["invalid_smiles"] += 1
            continue

        try:
            mol = Chem.MolFromSmiles(s, sanitize=True)
        except Exception:
            mol = None

        if mol is None:
            stats["invalid_smiles"] += 1
            continue

        if not passes_atom_set(mol):
            stats["atom_set_fail"] += 1
            continue

        if not is_single_fragment(mol):
            stats["fragment_fail"] += 1
            continue

        if not passes_size_limit(mol):
            stats["size_fail"] += 1
            continue

        if not passes_valence_rules(mol):
            stats["valence_fail"] += 1
            continue

        kept.append(Chem.MolToSmiles(mol, canonical=True))
        stats["kept_count"] += 1

    return kept, stats


def safe_get(scores, key):
    return scores[key] if key in scores else "N/A"


def main():
    RDLogger.DisableLog("rdApp.*")

    config = get_config()

    train_ds, eval_ds, test_ds, _ = datasets.get_dataset(config)
    train_smiles = train_ds.sub_smiles
    test_smiles = test_ds.sub_smiles
    test_scaffolds_smiles = test_smiles

    experiments = [
        {
            "name": "t25",
            "input_pkl": "exp/vpsde_qm9_cdgs/eval_t25_none_5k/pc_ckpt_200.pkl",
            "output_dir": "exp/vpsde_qm9_cdgs/eval_t25_filter_5k",
        },
        {
            "name": "t50",
            "input_pkl": "exp/vpsde_qm9_cdgs/eval_t50_none_5k/pc_ckpt_200.pkl",
            "output_dir": "exp/vpsde_qm9_cdgs/eval_t50_filter_5k",
        },
        {
            "name": "t100",
            "input_pkl": "exp/vpsde_qm9_cdgs/eval_t100_none_5k/pc_ckpt_200.pkl",
            "output_dir": "exp/vpsde_qm9_cdgs/eval_t100_filter_5k",
        },
        {
            "name": "t200",
            "input_pkl": "exp/vpsde_qm9_cdgs/eval_t200_none_5k/pc_ckpt_200.pkl",
            "output_dir": "exp/vpsde_qm9_cdgs/eval_t200_filter_5k",
        },
    ]

    for exp in experiments:
        name = exp["name"]
        input_pkl = exp["input_pkl"]
        output_dir = exp["output_dir"]

        if not os.path.exists(input_pkl):
            print(f"[{name}] missing input: {input_pkl}")
            continue

        os.makedirs(output_dir, exist_ok=True)

        with open(input_pkl, "rb") as f:
            raw_smiles = pickle.load(f)

        filtered_smiles, filter_stats = sanitize_and_filter(raw_smiles)
        unique_filtered_smiles = list(dict.fromkeys(filtered_smiles))

        train_set = set(train_smiles)
        novelty = (
            sum(1 for s in unique_filtered_smiles if s not in train_set) / len(unique_filtered_smiles)
            if len(unique_filtered_smiles) > 0 else 0.0
        )

        unique_ratio = (
            len(unique_filtered_smiles) / len(filtered_smiles)
            if len(filtered_smiles) > 0 else 0.0
        )

        metric_kwargs = {
            "gen": unique_filtered_smiles,
            "k": len(unique_filtered_smiles),
            "device": config.device,
            "n_jobs": 1,
            "test": test_smiles,
            "test_scaffolds": test_scaffolds_smiles,
            "train": train_smiles,
        }

        allowed = set(signature(get_all_metrics).parameters.keys())
        metric_kwargs = {k: v for k, v in metric_kwargs.items() if k in allowed}

        if len(unique_filtered_smiles) == 0:
            scores = {}
        else:
            try:
                scores = get_all_metrics(**metric_kwargs)
            except Exception as e:
                print(f"[{name}] get_all_metrics failed: {e}")
                scores = {}

        dropped_count = len(raw_smiles) - len(filtered_smiles)

        summary = {
            "input_count": len(raw_smiles),
            "filtered_count": len(filtered_smiles),
            "unique_filtered_count": len(unique_filtered_smiles),
            "dropped_count": dropped_count,
            "retention_ratio": len(filtered_smiles) / len(raw_smiles) if len(raw_smiles) > 0 else 0.0,
            "unique_ratio": unique_ratio,
            "valid": safe_get(scores, "valid"),
            "unique": safe_get(scores, f"unique@{len(unique_filtered_smiles)}"),
            "FCD/Test": safe_get(scores, "FCD/Test"),
            "Novelty_self": novelty,
            "invalid_smiles": filter_stats["invalid_smiles"],
            "atom_set_fail": filter_stats["atom_set_fail"],
            "fragment_fail": filter_stats["fragment_fail"],
            "size_fail": filter_stats["size_fail"],
            "valence_fail": filter_stats["valence_fail"],
        }

        with open(os.path.join(output_dir, "filtered_smiles.pkl"), "wb") as f:
            pickle.dump(unique_filtered_smiles, f)

        with open(os.path.join(output_dir, "summary.txt"), "w") as f:
            for k, v in summary.items():
                f.write(f"{k}: {v}\n")

        print(f"[{name}] done")
        for k, v in summary.items():
            print(f"  {k}: {v}")
        print()


if __name__ == "__main__":
    main()