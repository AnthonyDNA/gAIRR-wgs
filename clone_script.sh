# Find the site-packages directory of the activated Conda environment
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")
echo "Site-packages directory: $SITE_PACKAGES"

# The path to the package to be copied
SOURCE="gAIRR_suite"

rm -rf "$SITE_PACKAGES/gAIRR_suite"
cp -r "$SOURCE" "$SITE_PACKAGES/"