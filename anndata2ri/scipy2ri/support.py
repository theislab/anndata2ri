from typing import Union, Iterable, FrozenSet


supported_r_matrix_types = frozenset({"d", "l", "p"})
supported_r_matrix_storage = frozenset({"C", "R", "T", "di"})


def supported_r_matrix_classes(
    types: Union[Iterable[str], str] = supported_r_matrix_types,
    storage: Union[Iterable[str], str] = supported_r_matrix_storage,
) -> FrozenSet[str]:
    types = {types} if isinstance(types, str) else set(types)
    storage = {storage} if isinstance(storage, str) else set(storage)
    classes = {f"{t}g{s}Matrix" for t in types for s in storage - {"di"}}
    if "di" in storage:
        classes |= {f"{t}diMatrix" for t in types - {"p"}}
    return frozenset(classes)
