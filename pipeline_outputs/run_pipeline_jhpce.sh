#!/bin/bash
#$ -l bluejay,mem_free=60G,h_vmem=60G,h_fsize=800G
#$ -N SPEAQeasy_ex
#$ -o ./SPEAQeasy_ex.log
#$ -e ./SPEAQeasy_ex.log
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow /dcl01/lieber/ajaffe/Nick/RNAsp/main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "reverse" \
    --coverage \
    --annotation "/dcl01/lieber/ajaffe/Nick/RNAsp/Annotation" \
    -with-report execution_reports/JHPCE_run.html \
    -with-dag execution_DAGs/JHPCE_run.html \
    -profile jhpce \
    -w "/dcl01/lieber/ajaffe/lab/RNAsp_work/runs" \
    --input "/dcl01/lieber/ajaffe/lab/SPEAQeasy-example/sample_selection" \
    --output "/dcl01/lieber/ajaffe/lab/SPEAQeasy-example/pipeline_outputs"

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
#echo "Generating per-sample logs for debugging..."
#bash /dcl01/lieber/ajaffe/Nick/RNAsp/scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
