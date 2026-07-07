# Root Makefile for DiskMINT

# 1. Extract the version number from pyproject.toml
VERSION := $(shell grep -m 1 'version =' pyproject.toml | cut -d '"' -f 2)

# 2. Get the absolute path to your binaries
BIN_DIR := $(shell pwd)/chemistry/bin
REPO_DIR := $(shell pwd)

# .PHONY tells Make these are commands, not files
.PHONY: all chemistry install install_python clean uninstall uninstall_python link_bin check_env

# Default target: just compile Fortran
all: chemistry

# -------------------------------------------------------------
# Compilation & Installation Targets
# -------------------------------------------------------------

# 1. Compile Fortran (Calls the inner Makefile)
chemistry:
	@$(MAKE) -C chemistry/src
	@echo "----------------------------------------"
	@echo " ✅ Fortran Module Compiled Successfully"
	@echo "----------------------------------------"

# 2. Install Python Package (in editable mode)
install_python: check_env
	@python -m pip install -e .
	@echo "----------------------------------------"
	@echo " ✅ Python Module Installed Successfully"
	@echo "----------------------------------------"

# Report which environment 'make install' is about to install into, and warn
# (but don't fail) if none is active. Running 'make install' outside an
# activated env is the usual cause of "pip: command not found" reports (only
# 'pip3' exists on the base/system Python).
check_env:
	@if [ -n "$$CONDA_DEFAULT_ENV" ]; then \
		echo "📦 Installing into conda environment: $$CONDA_DEFAULT_ENV"; \
	elif [ -n "$$VIRTUAL_ENV" ]; then \
		echo "📦 Installing into virtualenv: $$VIRTUAL_ENV"; \
	else \
		echo "⚠️  No active conda environment or virtualenv detected."; \
		echo "   Run 'conda activate diskmint_stable' first (see environment.yml), then re-run 'make install'."; \
	fi
	@echo "   Python: $$(command -v python 2>/dev/null || command -v python3 2>/dev/null || echo 'not found')"

# 3. Main Install Command (Python + Fortran)
install: install_python chemistry
	@echo ""
	@echo "========================================"
	@echo " 🎉 DiskMINT v$(VERSION) Ready! 🎉"
	@echo "========================================"
	@echo ""
	@echo "Summary:"
	@if [ -n "$$CONDA_DEFAULT_ENV" ]; then \
		echo "  🐍 Python module   -> conda env '$$CONDA_DEFAULT_ENV', editable from $(REPO_DIR)"; \
	elif [ -n "$$VIRTUAL_ENV" ]; then \
		echo "  🐍 Python module   -> virtualenv '$$VIRTUAL_ENV', editable from $(REPO_DIR)"; \
	else \
		echo "  🐍 Python module   -> $$(command -v python 2>/dev/null || command -v python3 2>/dev/null), editable from $(REPO_DIR)"; \
	fi
	@echo "  ⚙️  Fortran binaries -> $(BIN_DIR)"
	@echo ""
	@echo "If you want to run the chemistry code from anywhere,"
	@echo ""
	@echo "👉 OPTION 1 (Recommended): Add to PATH"
	@echo "   Run this command (and add to ~/.zshrc if you use zsh):"
	@echo "   export PATH=\$$PATH:$(BIN_DIR)"
	@echo ""
	@echo "👉 OPTION 2: Create Shortcuts (Symlinks)"
	@echo "   Run 'make link_bin' to create shortcuts for Fortran-chemistry compiler in ~/.local/bin"
	@echo ""

# 4. Symlink Target: Mimics 'make install' of RADMC-3D
link_bin:
	@mkdir -p $(HOME)/.local/bin
	@ln -sf $(BIN_DIR)/disk_main $(HOME)/.local/bin/disk_main
	@ln -sf $(BIN_DIR)/disk_extract $(HOME)/.local/bin/disk_extract
	@echo "✅ Symlinks created in $(HOME)/.local/bin"
	@echo "   (Ensure ~/.local/bin is in your PATH)"

# -------------------------------------------------------------
# Cleanup & Uninstall Targets
# -------------------------------------------------------------

# Clean up compiled files
# (Also removes binaries from chemistry/bin to ensure a full clean)
clean:
	$(MAKE) -C chemistry/src clean
	@rm -f chemistry/bin/disk_main chemistry/bin/disk_extract
	@echo "Cleaned up compiled files and binaries."

# Uninstall Python only
uninstall_python:
	@python -m pip uninstall -y diskmint
	@echo "Uninstalled DiskMINT Python package."

# Full Uninstallation
# (Removes Python package, compiled files, and local symlinks)
uninstall: uninstall_python clean
	@rm -f $(HOME)/.local/bin/disk_main $(HOME)/.local/bin/disk_extract
	@echo "DiskMINT uninstalled (and symlinks removed)."