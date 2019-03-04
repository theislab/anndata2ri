dimension_reductions = {
    "pca",
    "dca",
    "tsne",  # There’s both tSNE and TSNE: https://rdrr.io/bioc/singleCellTK/man/getPCA.html
    "umap",
    ("diffmap", "dm"),  # X_diffmap is the official scanpy name
    "magic",
    "phate",
}
dimension_reductions_scanpy2sce = dict((n, n) if isinstance(n, str) else n for n in dimension_reductions)
dimension_reductions_sce2scanpy = {v: k for k, v in dimension_reductions_scanpy2sce.items()}


def scanpy2sce(name_scanpy: str) -> str:
    if name_scanpy.startswith("X_"):
        n = dimension_reductions_scanpy2sce.get(name_scanpy[2:])
        if n is not None:
            return n.upper()
    return name_scanpy


def sce2scanpy(name_sce: str) -> str:
    n = dimension_reductions_sce2scanpy.get(name_sce.lower())
    if n is not None:
        return f"X_{n}"
    return name_sce
