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
    return await MoleculeDAO.find_all_molecules(**request_body.to_dict())


@router.get(
        "/search",
        summary="Substructure search for molecules",
        response_description="List of molecules that match the substructure"
        )
async def search_molecule(substructure_smiles: str):
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
    check = await MoleculeDAO.update_molecule(mol_id=mol_id, **molecule.model_dump())
    if check:
        return {"message": f"The molecule with id {mol_id} has been updated!"}
    else:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Molecule Not Found"
            )
