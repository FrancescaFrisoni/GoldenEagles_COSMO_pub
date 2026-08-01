"""module with wrappers for dyntools

"""
from lagranto.tools import run_cmd


__all__ = ["run_cmd", "DynToolsException"]


class DynToolsException(RuntimeError):
    """Raise this when a dyn_tools call fails"""
