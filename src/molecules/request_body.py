class RBMolecule:
    def __init__(
        self,
        mol_id: int | None = None,
        name: str | None = None,
        smiles: str | None = None,
    ):
        self.mol_id = mol_id
        self.name = name
        self.smiles = smiles

    def to_dict(self) -> dict:
        data = {
            "mol_id": self.mol_id,
            "name": self.name,
            "smiles": self.smiles,
        }
        # Dict copy to avoid dict change while iteration
        filtered_data = {key: value for key, value in data.items() if value is not None}
        return filtered_data
