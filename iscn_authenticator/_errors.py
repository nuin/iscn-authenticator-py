# iscn_authenticator/_errors.py
"""Exception types shared by parser modules.

Re-exported from ``iscn_authenticator.parser`` for backwards compatibility.
"""


class ParseError(Exception):
    """Raised when a karyotype string cannot be parsed."""
