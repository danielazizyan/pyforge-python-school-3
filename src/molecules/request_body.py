from pydantic import BaseModel, Field
from typing import Optional


class RBMolecule(BaseModel):

    """
    Pydantic model for the request body when working with Molecules.

    Attributes:
        mol_id (Optional[int]): The molecule ID (optional).
        name (Optional[str]): The name of the molecule.
        smiles (Optional[str]): The SMILES string representation of the molecule.
    """

    mol_id: Optional[int] = Field(
                                None,
                                description="The molecule ID"
                            )
    name: Optional[str] = Field(
                                None,
                                description="The name of the molecule",
                                max_length=100
                            )
    smiles: Optional[str] = Field(
                                None,
                                description="The SMILES representation of the molecule",
                                max_length=255
                            )
