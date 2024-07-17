from rdkit import Chem

def substructure_search(mols, mol):
    """
    Search for molecules containing a given substructure.

    Parameters:
    mols (list of str): List of molecules as SMILES strings
    substructure (str): Substructure SMILES string

    Returns:
    list of str: List of molecule SMILES strings containing the given substructure.
    """

    matches = []
    substructure = Chem.MolFromSmiles(mol)

    if substructure is None:
        raise ValueError(f"Invalid SMILES string: {mol}")

    for smiles in mols:
        molecule = Chem.MolFromSmiles(smiles)

        if molecule is None:
            print(f"Invalid SMILES string: {smiles}")
            continue

        if molecule.HasSubstructMatch(substructure):
            matches.append(smiles)

    return matches

print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")) 
#["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]


