"""
Custom exception classes
"""


class AlreadyUploadedError(Exception):
    """Workflow and variant are already present in the database"""

    pass


class DatabaseCredentialError(Exception):
    """Missing or invalid database credentials"""

    pass


class VariantLookupError(Exception):
    """unrecoginised variant"""

    pass
