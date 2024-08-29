from __future__ import annotations

from functools import cache

from rpy2.robjects import Environment, packages


R_INT_BYTES = 4


@cache
def importr(name: str) -> packages.Package:
    return packages.importr(name)


@cache
def data(package: str, name: str | None = None) -> packages.PackageData | Environment:
    if name is None:
        return packages.data(importr(package))
    # Use cached version of PackageData collection and just fetch
    return data(package).fetch(name)
