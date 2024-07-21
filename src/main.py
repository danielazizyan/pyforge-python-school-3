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
            raise ValueError(f"Invalid SMILES string: {smiles}")

        if molecule.HasSubstructMatch(substructure):
            matches.append(smiles)

    return matches