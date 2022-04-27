class Error(Exception):
    """ limaps basic error"""

class IndividualNotProcessError(Error):
    """ Indivdual are not processed"""

class IndividualIntervalNotConsistent(Error):
    """ The interval of individuals are not consistent."""