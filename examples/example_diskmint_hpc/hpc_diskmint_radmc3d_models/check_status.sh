#!/bin/bash
# Monitor status of DiskMINT model runs

echo "==================================================="
echo "DiskMINT Model Run Status"
echo "==================================================="
echo ""

# Check SLURM queue
echo "Active SLURM jobs:"
squeue -u $USER --format="%.18i %.9P %.30j %.8T %.10M %.6D %R" | grep diskmint || echo "No active jobs"
echo ""

# Summary from logs
echo "Job completion summary:"
if [ -f job_submissions.log ]; then
    tail -n +2 job_submissions.log | while IFS=',' read -r mstar grid mark jobid nmodels; do
        echo "  Mstar=${mstar} (${grid}): JobID ${jobid}, ${nmodels} models"
        
        # Count completed/failed
        COMPLETED=$(find logs -name "diskmint_${jobid}_*.out" -type f 2>/dev/null | wc -l)
        FAILED=$(find logs -name "diskmint_${jobid}_*.err" -type f ! -size 0 2>/dev/null | wc -l)
        
        echo "    Completed: ${COMPLETED}/${nmodels}"
        if [ $FAILED -gt 0 ]; then
            echo "    Failed: ${FAILED}"
        fi
    done
else
    echo "No job submission log found"
fi
echo ""

# Check output directories
echo "Generated model outputs:"
if [ -d output ]; then
    N_OUTPUT=$(find output -maxdepth 1 -type d -name "diskmintmodel*" 2>/dev/null | wc -l)
    echo "  Total model directories: ${N_OUTPUT}"
    
    # Show recent completions
    echo "  Most recent (last 5):"
    find output -maxdepth 1 -type d -name "diskmintmodel*" -printf "%T+ %p\n" 2>/dev/null | \
        sort -r | head -5 | awk '{print "    " $2}' | sed 's|output/||'
else
    echo "  No output directory found"
fi
echo ""

# Check for errors in recent logs
echo "Recent errors (if any):"
RECENT_ERRORS=$(find logs -name "*.err" -type f ! -size 0 -mtime -1 2>/dev/null)
if [ -n "$RECENT_ERRORS" ]; then
    echo "$RECENT_ERRORS" | head -5 | while read errfile; do
        echo "  ${errfile}:"
        tail -3 "$errfile" | sed 's/^/    /'
    done
else
    echo "  None in last 24 hours"
fi
echo ""

echo "==================================================="
