from pydantic import BaseModel, Field, field_validator
from rdkit import Chem

SMILESField = Field(
    ..., min_length=1, max_length=100,
    description="SMILES representation of the molecule"
    )


class MoleculeBase(BaseModel):
    name: str = Field(..., min_length=1, max_length=100, description="Molecule Name")
    smiles: str = SMILESField

    @field_validator("smiles")
    def validate_smiles(cls, value):
        if Chem.MolFromSmiles(value) is None:
            raise ValueError("This SMILES has an invalid structure")
        return value

    model_config = {'from_attributes': True}


class MoleculeResponse(MoleculeBase):
    mol_id: int


class MoleculeAdd(MoleculeBase):
    pass


class MoleculeUpdate(MoleculeBase):
    pass
