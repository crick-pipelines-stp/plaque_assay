"""
Custom exception classes
"""


class AlreadyUploadedError(Exception):
    """Workflow and variant are already present in the database"""

    pass


class DatabaseCredentialError(Exception):
    """Missing or invalid database credentials"""

    pass


class VariantError(Exception):
    """unrecoginised variant"""

    pass
