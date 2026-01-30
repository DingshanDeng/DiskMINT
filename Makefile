# Root Makefile for DiskMINT

# 1. Extract the version number from pyproject.toml
#    (Greps the line 'version = "x.x.x"' and cuts out the quotes)
VERSION := $(shell grep -m 1 'version =' pyproject.toml | cut -d '"' -f 2)

# .PHONY tells Make these are commands, not files
.PHONY: all chemistry install install_python clean

# Default target: just compile Fortran
all: chemistry

# 2. Compile Fortran (Calls the inner Makefile)
chemistry:
	@$(MAKE) -C chemistry/src
	@echo ""
	@echo "----------------------------------------"
	@echo " ✅ Fortran Module Compiled Successfully"
	@echo "----------------------------------------"

# 3. Install Python Package (in editable mode)
install_python:
	@pip install -e .
	@echo ""
	@echo "----------------------------------------"
	@echo " ✅ Python Module Installed Successfully"
	@echo "----------------------------------------"

# 4. Full Installation (Python + Fortran)
install: install_python chemistry
	@echo ""
	@echo "========================================"
	@echo " 🎉 DiskMINT v$(VERSION) Ready for Use! 🎉"
	@echo "========================================"

# Clean up compiled files
clean:
	$(MAKE) -C chemistry/src clean
	@echo "Cleaned up compiled files."

# Uninstall
uninstall_python:
	@pip uninstall -y diskmint
	@echo "Uninstalled DiskMINT Python package."

# Full Uninstallation
uninstall: uninstall_python clean
	@echo "DiskMINT uninstalled."