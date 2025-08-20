#!/bin/bash
set -e  # stop if any command fails

# Render the site
quarto render

# Stage updated files
git add index.qmd docs/index.html

# Commit with timestamp message
git commit -m "Deploy site on $(date '+%Y-%m-%d %H:%M:%S')"

# Push to GitHub
git push origin main