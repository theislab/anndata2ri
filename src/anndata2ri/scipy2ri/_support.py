from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from collections.abc import Iterable


# these are documented in __init__.py because of sphinx limitations
supported_r_matrix_types = frozenset({'d', 'l', 'n'})
supported_r_matrix_storage = frozenset({'C', 'R', 'T', 'di'})


@lru_cache(maxsize=None)
def supported_r_matrix_classes(
    types: Iterable[str] | str = supported_r_matrix_types,
    storage: Iterable[str] | str = supported_r_matrix_storage,
) -> frozenset[str]:
    """Get supported classes, possibly limiting data types or storage types.

    :param types: Data type character(s) from :data:`supported_r_matrix_types`
    :param storage: Storage mode(s) from :data:`supported_r_matrix_storage`
    :return: All supported classes with those characters
    """
    types = {types} if isinstance(types, str) else set(types)
    storage = {storage} if isinstance(storage, str) else set(storage)

    bad_types = types - supported_r_matrix_types
    if bad_types:
        msg = f'Type(s) {bad_types} not supported.'
        raise ValueError(msg)
    bad_storage = storage - supported_r_matrix_storage
    if bad_storage:
        msg = f'Storage type(s) {bad_storage} not supported.'
        raise ValueError(msg)

    classes = {f'{t}g{s}Matrix' for t in types for s in storage - {'di'}}
    if 'di' in storage:
        classes |= {f'{t}diMatrix' for t in types - {'n'}}
    return frozenset(classes)
