"""
A set of wrappers for comman Matlab plotting commands so you don't have to use the m. prefix
"""

import pyspinw
m = pyspinw.Matlab()

def plot(*args, **kwargs):
    """
    Wrapper around Matlab plot() function
    """
    return m.plot(*args, **kwargs)

def subplot(*args, **kwargs):
    """
    Wrapper around Matlab subplot() function
    """
    return m.subplot(*args, **kwargs)

def xlim(*args, **kwargs):
    """
    Wrapper around Matlab xlim() function
    """
    if args and isinstance(args[0], str):
        return m.xlim(*args, **kwargs)
    else:
        m.xlim(*args, **kwargs)

def ylim(*args, **kwargs):
    """
    Wrapper around Matlab ylim() function
    """
    if args and isinstance(args[0], str):
        return m.ylim(*args, **kwargs)
    else:
        m.ylim(*args, **kwargs)

def xlabel(*args, **kwargs):
    """
    Wrapper around Matlab xlabel() function
    """
    return m.xlabel(*args, **kwargs)

def ylabel(*args, **kwargs):
    """
    Wrapper around Matlab ylabel() function
    """
    return m.ylabel(*args, **kwargs)

def set(*args, **kwargs):
    """
    Wrapper around Matlab set() function
    """
    m.set(*args, **kwargs)

def get(*args, **kwargs):
    """
    Wrapper around Matlab get() function
    """
    return m.get(*args, **kwargs)

def gca(*args, **kwargs):
    """
    Wrapper around Matlab gca() function
    """
    return m.gca(*args, **kwargs)

def gcf(*args, **kwargs):
    """
    Wrapper around Matlab gcf() function
    """
    return m.gcf(*args, **kwargs)

def legend(*args, **kwargs):
    """
    Wrapper around Matlab legend() function
    """
    return m.legend(*args, **kwargs)

def hold(*args, **kwargs):
    """
    Wrapper around Matlab hold() function
    """
    m.hold(*args, **kwargs)

def pcolor(*args, **kwargs):
    """
    Wrapper around Matlab pcolor() function
    """
    return m.pcolor(*args, **kwargs)

def contour(*args, **kwargs):
    """
    Wrapper around Matlab contour() function
    """
    return m.contour(*args, **kwargs)

