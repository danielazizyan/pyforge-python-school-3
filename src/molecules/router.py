from fastapi import APIRouter, Depends, HTTPException, status
from src.molecules.dao import MoleculeDAO
from src.molecules.request_body import RBMolecule
from src.molecules.schema import MoleculeResponse, MoleculeAdd, MoleculeUpdate
from rdkit import Chem

router = APIRouter(prefix="/molecules", tags=["molecules"])


@router.get(
        "/",
        summary="Get all molecules",
        response_description="List of all molecules"
        )
async def get_all_molecules(
    request_body: RBMolecule = Depends(),
) -> list[MoleculeResponse]:

    """
    Retrieve a list of all molecules in the database.

    Parameters:
    - request_body: A request body model (RBMolecule) containing filters for the query.

    Returns:
    - A list of MoleculeResponse objects that match the query filters.
    """

    return await MoleculeDAO.find_all_molecules(**request_body.to_dict())


@router.get(
        "/search",
        summary="Substructure search for molecules",
        response_description="List of molecules that match the substructure"
        )
async def search_molecule(substructure_smiles: str):

    """
    Search for molecules containing a specific substructure.

    Parameters:
    - substructure_smiles: The SMILES string representing the substructure to search for.

    Returns:
    - A list of molecules that contain the given substructure.

    Raises:
    - HTTPException: If the SMILES string is invalid or if no molecules are found.
    """

    try:
        if Chem.MolFromSmiles(substructure_smiles) is None:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="This SMILES has an invalid structure"
                )

        matches = await MoleculeDAO.search_molecule(
            substructure_smiles=substructure_smiles
            )
        if matches:
            return matches
        else:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="No molecules found with the given substructure"
            )
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.get(
        "/{mol_id}",
        summary="Get a molecule by ID",
        response_description="The molecule data"
        )
async def get_molecule_by_id(mol_id: int) -> MoleculeResponse | dict:

    """
    Retrieve a molecule by its ID.

    Parameters:
    - mol_id: The ID of the molecule to retrieve.

    Returns:
    - A MoleculeResponse object containing the molecule's data, or a dictionary.

    Raises:
    - HTTPException: If the molecule with the specified ID is not found.
    """

    rez = await MoleculeDAO.find_full_data(mol_id=mol_id)
    if rez is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Molecule not found"
            )
    return rez


@router.post(
        "/add/",
        summary="Add a new molecule",
        response_description="The added molecule data"
        )
async def add_molecule(molecule: MoleculeAdd) -> dict:

    """
    Add a new molecule to the database.

    Parameters:
    - molecule: A MoleculeAdd object containing the data for the new molecule.

    Returns:
    - A confirmation message and the data of the added molecule.

    Raises:
    - HTTPException: If there is an error adding the molecule.
    """

    check = await MoleculeDAO.add_molecule(**molecule.model_dump())
    if check:
        return {"message": "The molecule is added!", "molecule": molecule}
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Error adding a molecule"
            )


@router.delete(
        "/delete/{mol_id}",
        summary="Delete a molecule by ID",
        response_description="Confirmation of deletion"
        )
async def delete_molecule_by_id(mol_id: int) -> dict:

    """
    Delete a molecule from the database by its ID.

    Parameters:
    - mol_id: The ID of the molecule to delete.

    Returns:
    - A confirmation message indicating the molecule has been deleted.

    Raises:
    - HTTPException: If the molecule with the specified ID is not found.
    """

    check = await MoleculeDAO.delete_molecule_by_id(mol_id=mol_id)
    if check:
        return {"message": f"The molecule with id {mol_id} is deleted!"}
    else:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Molecule Not Found"
            )


@router.put(
        "/{mol_id}",
        summary="Update a molecule by ID",
        response_description="Confirmation of update"
        )
async def update_molecule(mol_id: int, molecule: MoleculeUpdate) -> dict:
    """
    Update the data of a molecule by its ID.

    Parameters:
    - mol_id: The ID of the molecule to update.
    - molecule: A MoleculeUpdate object containing the updated molecule data.

    Returns:
    - A confirmation message indicating the molecule has been updated.

    Raises:
    - HTTPException: If the molecule with the specified ID is not found.
    """
    check = await MoleculeDAO.update_molecule(mol_id=mol_id, **molecule.model_dump())
    if check:
        return {"message": f"The molecule with id {mol_id} has been updated!"}
    else:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Molecule Not Found"
            )
