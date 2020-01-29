from functools import lru_cache
from typing import Optional, Union

from rpy2.robjects import packages, Environment


@lru_cache()
def importr(name: str) -> packages.Package:
    return packages.importr(name)


@lru_cache()
def data(package: str, name: Optional[str] = None) -> Union[packages.PackageData, Environment]:
    if name is None:
        return packages.data(importr(package))
    else:
        # Use cached version of PackageData collection and just fetch
        return data(package).fetch(name)
