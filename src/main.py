from os import getenv
from fastapi import FastAPI, status, HTTPException
from rdkit import Chem

app = FastAPI()

smiles_db = [
    {"mol_id": 1, "name": "Ethanol", "smiles": "CCO"},
    {"mol_id": 2, "name": "Benzene", "smiles": "c1ccccc1"},
    {"mol_id": 3, "name": "Acetic acid", "smiles": "CC(=O)O"},
    {"mol_id": 4, "name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
    {"mol_id": 5, "name": "Methane", "smiles": "C"},
    {"mol_id": 6, "name": "Propane", "smiles": "CCC"},
    {"mol_id": 7, "name": "Toluene", "smiles": "Cc1ccccc1"},
    {"mol_id": 8, "name": "Phenol", "smiles": "c1ccc(cc1)O"},
    {"mol_id": 9, "name": "Formaldehyde", "smiles": "C=O"},
    {"mol_id": 10, "name": "Acetone", "smiles": "CC(=O)C"},
    {"mol_id": 11, "name": "Cyclohexane", "smiles": "C1CCCCC1"},
    {"mol_id": 12, "name": "Ethylene", "smiles": "C=C"},
    {"mol_id": 13, "name": "Styrene", "smiles": "c1ccccc1C=C"},
    {"mol_id": 14, "name": "Hexane", "smiles": "CCCCCC"}
]


def validate_smiles(smiles: str):

    """
    Validates a SMILES string to ensure it can be correctly parsed into a molecule.

    Parameters:
    - smiles: The SMILES string to be validated.

    Raises:
    - ValueError: If the SMILES string cannot be parsed into a molecule.
    """

    if not smiles or Chem.MolFromSmiles(smiles) is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")


@app.get("/")
def get_server():

    """
    Retrieve the server ID.
    """

    return {"server_id": getenv("SERVER_ID", "1")}


@app.get(
        "/molecules",
        tags=["Molecules"],
        summary="List all molecules",
        response_description="A list of all molecules in the database"
        )
def list_molecules():

    """
    List all molecules.
    """

    return smiles_db


@app.post(
        "/add",
        status_code=status.HTTP_201_CREATED,
        tags=["Molecules"], summary="Add a new molecule",
        response_description="Molecule has been added"
        )
def add_molecule(molecule: dict):

    """
    Add molecule (smiles) to the database.

    Parameters:
    - molecule: A dictionary containing mol_id, name, and smiles.

    Returns:
    - The created molecule.
    """

    for mol in smiles_db:
        if mol["mol_id"] == molecule["mol_id"]:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Molecule with this ID already exists."
                )

    try:
        validate_smiles(molecule["smiles"])
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))

    smiles_db.append(molecule)
    return molecule


@app.get(
        "/molecules/{mol_id}",
        tags=["Molecules"],
        summary="Get a molecule by ID",
        response_description="The molecule has been retrieved"
        )
def get_molecule(mol_id: int):

    """
    Get molecule by identifier.

    Parameters:
    - mol_id: The ID of the molecule to retrieve.

    Returns:
    - The molecule with the specified ID.
    """

    for molecule in smiles_db:
        if molecule["mol_id"] == mol_id:
            return molecule
    raise HTTPException(
        status_code=status.HTTP_404_NOT_FOUND, detail="Molecule Not Found"
        )


@app.put(
        "/molecules/{mol_id}",
        tags=["Molecules"],
        summary="Update a molecule by ID",
        response_description="The molecule has been updated"
        )
def update_molecule(mol_id: int, updated_molecule: dict):

    """
    Updating a molecule by identifier.

    Parameters:
    - mol_id: The ID of the molecule to update.
    - updated_molecule: A dictionary containing the updated molecule data.

    Returns:
    - The updated molecule.
    """

    for index, molecule in enumerate(smiles_db):
        if molecule["mol_id"] == mol_id:

            try:
                validate_smiles(updated_molecule["smiles"])
            except ValueError as e:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST, detail=str(e)
                    )

            smiles_db[index] = updated_molecule
            return updated_molecule
    raise HTTPException(
        status_code=status.HTTP_404_NOT_FOUND, detail="Molecule Not Found"
        )


@app.delete(
        "/molecules/{mol_id}",
        tags=["Molecules"],
        summary="Delete a molecule by ID",
        response_description="Molecule has been deleted"
        )
def delete_molecule(mol_id: int):

    """
    Delete a molecule by identifier.

    Parameters:
    - mol_id: The ID of the molecule to delete.

    Returns:
    - The deleted molecule.
    """

    for index, molecule in enumerate(smiles_db):
        if molecule["mol_id"] == mol_id:
            deleted_molecule = smiles_db.pop(index)
            return deleted_molecule
    raise HTTPException(
        status_code=status.HTTP_404_NOT_FOUND, detail="Molecule Not Found"
        )


@app.get(
        "/search",
        tags=["Molecules"],
        summary="Substructure search for molecules",
        response_description=(
            "Retrieved a list of molecules containing the given substructure"
        )
        )
def search_molecule(substructure_smiles: str):

    """
    Substructure search for all added molecules.

    Parameters:
    - substructure_smiles: The substructure SMILES string to search for.

    Returns:
    - A list of molecules containing the given substructure.
    """

    try:
        smiles_list = [mol["smiles"] for mol in smiles_db]
        matching_smiles = substructure_search(smiles_list, substructure_smiles)
        matches = [mol for mol in smiles_db if mol["smiles"] in matching_smiles]
        return matches
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


def substructure_search(mols, substructure_smiles):
    """
    Search for molecules containing a given substructure.

    Parameters:
    mols (list of str): List of molecules as SMILES strings
    substructure (str): Substructure SMILES string

    Returns:
    list of str: List of molecule SMILES strings containing the given substructure.
    """

    validate_smiles(substructure_smiles)
    substructure = Chem.MolFromSmiles(substructure_smiles)

    matches = []
    for smiles in mols:
        validate_smiles(smiles)
        molecule = Chem.MolFromSmiles(smiles)

        if molecule.HasSubstructMatch(substructure):
            matches.append(smiles)

    return matches
