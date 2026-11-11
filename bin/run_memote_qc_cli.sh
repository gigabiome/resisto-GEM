#!/bin/bash
# Run memote quality control on all models using CLI
# This approach is often more stable than Python API

set -e

RESULTS_DIR="${1:-.}"
MODELS_DIR="$RESULTS_DIR/results/models"
QC_DIR="$RESULTS_DIR/results/memote_qc"
SUMMARY_FILE="$RESULTS_DIR/results/memote_summary.csv"

# Check memote is available
if ! command -v memote &> /dev/null; then
    echo "ERROR: memote not found. Install with: pip install memote"
    exit 1
fi

# Validate paths
if [ ! -d "$MODELS_DIR" ]; then
    echo "ERROR: Models directory not found: $MODELS_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$QC_DIR"

echo "========================================================================"
echo "Running memote quality control on all models"
echo "========================================================================"
echo ""

# Find all SBML files
MODEL_COUNT=$(find "$MODELS_DIR" -name "*.xml" -type f | wc -l)
echo "Found $MODEL_COUNT models to test"
echo ""

# Create summary header
echo "model,score,passed,warnings,failed,reactions,metabolites,genes,biomass_reactions" > "$SUMMARY_FILE"

# Counter
CURRENT=0

# Run memote on each model
for MODEL_FILE in "$MODELS_DIR"/*.xml; do
    CURRENT=$((CURRENT + 1))
    MODEL_NAME=$(basename "$MODEL_FILE" .xml)
    JSON_OUTPUT="$QC_DIR/$MODEL_NAME.json"

    printf "[%3d/%d] %-40s " "$CURRENT" "$MODEL_COUNT" "$MODEL_NAME"

    # Run memote and capture JSON output
    if memote run --model "$MODEL_FILE" --json "$JSON_OUTPUT" > /dev/null 2>&1; then
        # Extract metrics from JSON
        if [ -f "$JSON_OUTPUT" ]; then
            SCORE=$(python3 -c "import json; d=json.load(open('$JSON_OUTPUT')); print(f\"{d.get('score', 0):.3f}\")" 2>/dev/null || echo "0.000")

            # Count passed/failed tests (simplified)
            PASSED=$(python3 -c "
import json
d = json.load(open('$JSON_OUTPUT'))
tests = d.get('tests', {})
print(sum(1 for t in tests.values() if t.get('metric') == 1.0))
" 2>/dev/null || echo "0")

            FAILED=$(python3 -c "
import json
d = json.load(open('$JSON_OUTPUT'))
tests = d.get('tests', {})
print(sum(1 for t in tests.values() if t.get('metric') == 0.0))
" 2>/dev/null || echo "0")

            WARNINGS=$(python3 -c "
import json
d = json.load(open('$JSON_OUTPUT'))
tests = d.get('tests', {})
print(sum(1 for t in tests.values() if 0 < t.get('metric', 0) < 1.0))
" 2>/dev/null || echo "0")

            echo "✓ (${SCORE})"

            # Try to get model stats (using COBRApy)
            STATS=$(python3 << 'PYSCRIPT' 2>/dev/null || echo "0,0,0,0")
try:
    from cobra.io import read_sbml_model
    model = read_sbml_model("$MODEL_FILE")
    reactions = len(model.reactions)
    metabolites = len(model.metabolites)
    genes = len(model.genes)
    biomass = len([r for r in model.reactions if 'biomass' in r.id.lower()])
    print(f"{reactions},{metabolites},{genes},{biomass}")
except:
    print("0,0,0,0")
PYSCRIPT
            )

            # Append to summary
            echo "$MODEL_NAME,$SCORE,$PASSED,$WARNINGS,$FAILED,$STATS" >> "$SUMMARY_FILE"
        else
            echo "✗ (no output)"
        fi
    else
        echo "✗ (test failed)"
    fi
done

echo ""
echo "========================================================================"
echo "✓ Memote QC complete"
echo "========================================================================"
echo ""
echo "Results:"
echo "  JSON reports: $QC_DIR/"
echo "  Summary CSV:  $SUMMARY_FILE"
echo ""

# Print summary statistics
python3 << 'PYSCRIPT'
import csv

with open("$SUMMARY_FILE") as f:
    reader = csv.DictReader(f)
    scores = []
    passed = 0
    warning = 0
    failed = 0

    for row in reader:
        try:
            score = float(row['score'])
            scores.append(score)
            if score >= 0.90:
                passed += 1
            elif score >= 0.70:
                warning += 1
            else:
                failed += 1
        except:
            pass

    if scores:
        print(f"Quality summary:")
        print(f"  ✓ High quality (≥90%):   {passed}")
        print(f"  ⚠ Needs review (70-90%): {warning}")
        print(f"  ✗ Poor quality (<70%):   {failed}")
        print(f"\nScore statistics:")
        print(f"  Average: {sum(scores)/len(scores):.1%}")
        print(f"  Min:     {min(scores):.1%}")
        print(f"  Max:     {max(scores):.1%}")
PYSCRIPT
