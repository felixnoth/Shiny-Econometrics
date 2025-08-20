#!/bin/bash
set -euo pipefail

# Render site
quarto render

# Stage common site changes (+ apps and styles)
git add -A docs index.qmd _quarto.yml styles.css apps || true

# Commit only if there are staged changes
if git diff --cached --quiet; then
  echo "No changes to commit."
else
  git commit -m "Deploy site on $(date '+%Y-%m-%d %H:%M:%S %Z')"
  git push origin main
fi

# Cache-busted URL for quick check
echo "Open: https://felixnoth.github.io/Shiny-Econometrics/?v=$(date +%s)"