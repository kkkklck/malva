# Copyright (c) 2025 Daniel León-Periñán and Nikolaos Karaiskos
#                    Rajewsky Lab, Max Delbrück Center for Molecular Medicine (MDC), Berlin
#
# Non-commercial and academic use only. See LICENSE for full terms.

try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("malva")
except PackageNotFoundError:
    # Package not yet installed (e.g., running from source tree)
    __version__ = "unknown"

from malva.malva_index import MalvaIndex

__all__ = ["MalvaIndex"]
