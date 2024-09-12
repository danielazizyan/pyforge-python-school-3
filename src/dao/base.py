import logging
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.future import select
from sqlalchemy import update as sqlalchemy_update, delete as sqlalchemy_delete
from src.database import async_session_maker
from sqlalchemy.ext.asyncio import AsyncSession


logger = logging.getLogger(__name__)


class BaseDAO:

    """
    A generic Data Access Object (DAO) class providing common database operations
    such as adding, updating, deleting, and retrieving records.

    Attributes:
        model (SQLAlchemy model): The model associated with the DAO.
    """

    model = None

    @classmethod
    async def _commit(cls, session: AsyncSession):

        """
        Commits the current transaction. If an error occurs, rolls back the transaction.

        Parameters:
            session (AsyncSession): The current database session.

        Raises:
            SQLAlchemyError: If there is an error during the commit.
        """

        try:
            await session.commit()
        except SQLAlchemyError as e:
            await session.rollback()
            raise e

    @classmethod
    async def find_one_or_none_by_id(cls, data_id: int):

        """
        Finds a single record by its primary key ID. Returns None if no record is found.

        Parameters:
            data_id (int): The primary key ID of the record to retrieve.

        Returns:
            The record with the given ID, or None if no such record exists.
        """

        async with async_session_maker() as session:
            query = select(cls.model).filter_by(id=data_id)
            logger.debug(f"Executing query: {query}")
            result = await session.execute(query)
            return result.scalar_one_or_none()

    @classmethod
    async def find_one_or_none(cls, **filter_by):

        """
        Finds a single record based on the specified filters.
        Returns None if no record is found.

        Parameters:
            **filter_by: Keyword arguments to filter records by.

        Returns:
            The record that matches the filters, or None if no such record exists.
        """

        async with async_session_maker() as session:
            query = select(cls.model).filter_by(**filter_by)
            logger.debug(f"Executing query: {query}")
            result = await session.execute(query)
            return result.scalar_one_or_none()

    @classmethod
    async def find_all(cls, **filter_by):

        """
        Finds all records that match the given filters.

        Parameters:
            **filter_by: Keyword arguments to filter records by.

        Returns:
            A list of all records that match the filters.
        """

        async with async_session_maker() as session:
            query = select(cls.model).filter_by(**filter_by)
            logger.debug(f"Executing query: {query}")
            result = await session.execute(query)
            return result.scalars().all()

    @classmethod
    async def add(cls, **values):

        """
        Adds a new record to the database.

        Parameters:
            **values: The values to set on the new record.

        Returns:
            The newly created record.
        """

        async with async_session_maker() as session:
            async with session.begin():
                new_instance = cls.model(**values)
                session.add(new_instance)
                logger.debug(f"Adding instance: {new_instance}")
                await cls._commit(session)
                return new_instance

    @classmethod
    async def add_many(cls, instances: list[dict]):

        """
        Adds multiple new records to the database.

        Parameters:
            instances: A list of dictionaries representing the values for each record.

        Returns:
            A list of the newly created records.
        """

        async with async_session_maker() as session:
            async with session.begin():
                new_instances = [cls.model(**values) for values in instances]
                session.add_all(new_instances)
                logger.debug(f"Adding instances: {new_instances}")
                await cls._commit(session)
                return new_instances

    @classmethod
    async def update(cls, filter_by, **values):

        """
        Updates records that match the given filters with the provided values.

        Parameters:
            filter_by (dict): A dictionary of filters to select the records to update.
            **values: The values to update the records with.

        Returns:
            The number of records that were updated.
        """

        async with async_session_maker() as session:
            async with session.begin():
                query = (
                    sqlalchemy_update(cls.model)
                    .where(*[getattr(cls.model, k) == v for k, v in filter_by.items()])
                    .values(**values)
                    .execution_options(synchronize_session="fetch")
                )
                logger.debug(f"Executing query: {query}")
                result = await session.execute(query)
                await cls._commit(session)
                return result.rowcount

    @classmethod
    async def delete(cls, delete_all: bool = False, **filter_by):

        """
        Deletes records that match the given filters.
        If no filters are provided, raises an error.

        Parameters:
            delete_all (bool): If True, all records will be deleted.
            **filter_by: Keyword arguments to filter records by.

        Returns:
            The number of records that were deleted.

        Raises:
            ValueError: If `delete_all` is False and no filters are provided.
        """

        if delete_all is False:
            if not filter_by:
                raise ValueError(
                    "You need to provide at least one parameter to delete."
                )

        async with async_session_maker() as session:
            async with session.begin():
                query = sqlalchemy_delete(cls.model).filter_by(**filter_by)
                logger.debug(f"Executing query: {query}")
                result = await session.execute(query)
                await cls._commit(session)
                return result.rowcount
