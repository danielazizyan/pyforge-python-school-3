from sqlalchemy.orm import Mapped
from src.database import Base, int_pk, str_uniq


class Molecule(Base):
    mol_id: Mapped[int_pk]
    name: Mapped[str]
    smiles: Mapped[str_uniq]

    def __str__(self):
        return (
            f"{self.__class__.__name__}(mol_id={self.mol_id}, "
            f"name={self.name!r},"
            f"smiles={self.smiles!r})"
        )

    def __repr__(self):
        return str(self)
