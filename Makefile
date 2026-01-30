# Root Makefile for DiskMINT

# 1. Extract the version number from pyproject.toml
VERSION := $(shell grep -m 1 'version =' pyproject.toml | cut -d '"' -f 2)

# 2. Get the absolute path to your binaries
BIN_DIR := $(shell pwd)/chemistry/bin

# .PHONY tells Make these are commands, not files
.PHONY: all chemistry install install_python clean uninstall uninstall_python link_bin

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
install_python:
	@pip install -e .
	@echo "----------------------------------------"
	@echo " ✅ Python Module Installed Successfully"
	@echo "----------------------------------------"

# 3. Main Install Command (Python + Fortran)
install: install_python chemistry
	@echo ""
	@echo "========================================"
	@echo " 🎉 DiskMINT v$(VERSION) Ready! 🎉"
	@echo "========================================"
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
	@pip uninstall -y diskmint
	@echo "Uninstalled DiskMINT Python package."

# Full Uninstallation
# (Removes Python package, compiled files, and local symlinks)
uninstall: uninstall_python clean
	@rm -f $(HOME)/.local/bin/disk_main $(HOME)/.local/bin/disk_extract
	@echo "DiskMINT uninstalled (and symlinks removed)."