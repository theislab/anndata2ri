SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = AnnData2RI
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help Makefile

help: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
