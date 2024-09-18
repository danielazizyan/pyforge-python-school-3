from fastapi import APIRouter, Depends, HTTPException, status, Query
from src.molecules.dao import MoleculeDAO
from src.molecules.request_body import RBMolecule
from src.molecules.schema import MoleculeResponse, MoleculeAdd, MoleculeUpdate
from rdkit import Chem
import logging

router = APIRouter(prefix="/molecules", tags=["molecules"])

logger = logging.getLogger(__name__)


@router.get(
    "/",
    summary="Get all molecules",
    response_description="List of all molecules",
)
async def get_all_molecules(
    request_body: RBMolecule = Depends(),
    limit: int = Query(None, description="Limit the number of molecules returned")
) -> list[MoleculeResponse]:

    """
    Retrieve a list of all molecules in the database, with an optional limit.

    Parameters:
    - request_body: A request body model (RBMolecule) containing filters for the query.
    - limit: Optional parameter to limit the number of molecules returned.

    Returns:
    - A list of MoleculeResponse objects that match the query filters.
    """

    logger.info(f"Fetching all molecules with a limit of {limit}")
    molecule_list = []

    async for molecule in MoleculeDAO.find_all_molecules(limit=limit):
        molecule_list.append(MoleculeResponse(**molecule.__dict__))

    logger.info(f"Fetched {len(molecule_list)} molecules.")
    return molecule_list


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

    logger.info(f"Searching for molecules with substructure: {substructure_smiles}")
    try:
        if Chem.MolFromSmiles(substructure_smiles) is None:
            logger.warning(f"Invalid SMILES structure: {substructure_smiles}")
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="This SMILES has an invalid structure"
                )

        matches = await MoleculeDAO.search_molecule(
            substructure_smiles=substructure_smiles
            )
        if matches:
            logger.info(
                f"Found {len(matches)} molecules matching substructure"
                f"{substructure_smiles}."
            )
            return matches
        else:
            logger.info(f"No molecules found with substructure {substructure_smiles}.")
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="No molecules found with the given substructure"
            )
    except ValueError as e:
        logger.error(f"ValueError during substructure search: {e}")
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

    logger.info(f"Fetching molecule with ID: {mol_id}")
    rez = await MoleculeDAO.find_full_data(mol_id=mol_id)
    if rez is None:
        logger.warning(f"Molecule with ID {mol_id} not found.")
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Molecule not found"
            )
    logger.info(f"Molecule with ID {mol_id} found.")
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

    logger.info(f"Adding new molecule: {molecule.name}")
    check = await MoleculeDAO.add_molecule(**molecule.model_dump())
    if check:
        logger.info(f"Molecule {molecule.name} added successfully.")
        return {"message": "The molecule is added!", "molecule": molecule}
    else:
        logger.error(f"Error adding molecule {molecule.name}")
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

    logger.info(f"Deleting molecule with ID: {mol_id}")
    check = await MoleculeDAO.delete_molecule_by_id(mol_id=mol_id)
    if check:
        logger.info(f"Molecule with ID {mol_id} deleted successfully.")
        return {"message": f"The molecule with id {mol_id} is deleted!"}
    else:
        logger.warning(f"Molecule with ID {mol_id} not found for deletion.")
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

    logger.info(f"Updating molecule with ID: {mol_id}")
    check = await MoleculeDAO.update_molecule(mol_id=mol_id, **molecule.model_dump())
    if check:
        logger.info(f"Molecule with ID {mol_id} updated successfully.")
        return {"message": f"The molecule with id {mol_id} has been updated!"}
    else:
        logger.warning(f"Molecule with ID {mol_id} not found for update.")
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Molecule Not Found"
            )
