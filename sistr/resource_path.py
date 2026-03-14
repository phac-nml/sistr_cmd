"""Resource path helper using importlib.resources (replaces deprecated pkg_resources)."""
from __future__ import print_function

from importlib.resources import files


def resource_path(package: str, path: str) -> str:
    """Get filesystem path to a package resource (replaces pkg_resources.resource_filename)."""
    parts = path.replace('\\', '/').strip('/').split('/')
    result = files(package)
    for part in parts:
        if part:
            result = result / part
    return str(result)
