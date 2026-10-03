SPHINXOPTS    ?=
SOURCEDIR     = .
BUILDDIR      = _build
VENV_DIR      = .venv
VENV_ACTIVATE = $(VENV_DIR)/bin/activate
PYTHON        = $(VENV_DIR)/bin/python3
SPHINXPROJ    = ORCA Python Interface Chinese documentation
SPHINXBUILD   = $(PYTHON) -m sphinx
SPHINXAUTO    = $(PYTHON) -m sphinx_autobuild

SHELL := /bin/bash

clean_api:
	rm -rf "$(SOURCEDIR)/contents/api"
	rm -rf $(SOURCEDIR)/contents/scripts/*.sha256

clean:
	rm -rf "$(BUILDDIR)/html"

distclean: clean clean_api
	rm -rf "$(BUILDDIR)"
	rm -rf $(VENV_DIR)/

$(VENV_DIR)/bin/uv:
	@if [ ! -d "$(VENV_DIR)/bin/activate" ]; then \
		python3 -m venv $(VENV_DIR); \
		source $(VENV_ACTIVATE) && $(PYTHON) -m pip install --upgrade pip; \
	fi
	@if [ ! -f "$(VENV_DIR)/bin/uv" ]; then \
		source $(VENV_ACTIVATE) && $(PYTHON) -m pip install --upgrade uv; \
	fi

sphinx: $(VENV_DIR)/bin/uv
	source $(VENV_ACTIVATE) && $(PYTHON) -m uv sync --active --group docs --no-editable --locked --no-python-downloads

auto: sphinx
	source $(VENV_ACTIVATE) && $(SPHINXAUTO) --re-ignore $(BUILDDIR) $(SOURCEDIR) $(BUILDDIR)/html

html: sphinx
	source $(VENV_ACTIVATE) && $(SPHINXBUILD) $(SPHINXOPTS) -b html . $(BUILDDIR)/html

# > Phony targets
.PHONY: auto clean clean_generated distclean html sphinx

# Additional output formats use the same translated source tree.
singlehtml epub latex: sphinx
	source $(VENV_ACTIVATE) && $(SPHINXBUILD) $(SPHINXOPTS) -b $@ . $(BUILDDIR)/$@

.PHONY: singlehtml epub latex
