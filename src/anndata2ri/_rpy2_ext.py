from __future__ import annotations

from functools import lru_cache

from rpy2.robjects import Environment, packages


@lru_cache
def importr(name: str) -> packages.Package:
    return packages.importr(name)


@lru_cache
def data(package: str, name: str | None = None) -> packages.PackageData | Environment:
    if name is None:
        return packages.data(importr(package))
    # Use cached version of PackageData collection and just fetch
    return data(package).fetch(name)
