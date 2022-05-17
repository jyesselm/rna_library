"""
Error classes to be inherited from.
"""


class InvalidDotBracket(Exception):
    """
	Exception for a mal-formed dot-bracket secondary structure.
	"""

    pass


class MissingDependency(Exception):
    """
    Exception for when a dependency is missing in the system.
    """

    pass


class InvalidArgument(Exception):
   """
   Exception for a bad argument being supplied.
   """
   pass
