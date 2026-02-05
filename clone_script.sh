# Find the site-packages directory of the activated Conda environment
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")
echo "Site-packages directory: $SITE_PACKAGES"

# The path to the package to be copied
SOURCE="gAIRR_suite"

rm -rf "$SITE_PACKAGES/gAIRR_suite"
cp -r "$SOURCE" "$SITE_PACKAGES/"

ENV_BIN="$(python -c 'import sys, os; print(os.path.dirname(sys.executable))')"
echo "Environment bin directory: $ENV_BIN"

cat > "$ENV_BIN/gAIRR_wgs" << 'PY'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re, sys
from gAIRR_suite.gAIRR_wgs import main

if __name__ == "__main__":
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    raise SystemExit(main())
PY

chmod +x "$ENV_BIN/gAIRR_wgs"


# Install gAIRR_extract
cat > "$ENV_BIN/gAIRR_extract" << 'SHELL'
#!/bin/bash
# Wrapper script to execute the gAIRR extract shell script from the installed package

SCRIPT_DIR="$(python -c "import site; print(site.getsitepackages()[0])")/gAIRR_suite"
exec bash "$SCRIPT_DIR/gAIRR_extract.sh" "$@"
SHELL

chmod +x "$ENV_BIN/gAIRR_extract"


