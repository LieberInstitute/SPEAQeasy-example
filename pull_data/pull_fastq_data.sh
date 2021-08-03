#  Run this script from the main 'SPEAQeasy-example' folder, specifying an
#  output folder to place downloaded files:
#
#  bash pull_data/pull_fastq_data.sh "/some_directory"

if [ "$1" == "" ]; then
    echo "Please provide a download location for the FASTQ data when running this script."
    echo -e '\ne.g. bash pull_data/pull_fastq_data.sh "/some_directory"'
    exit 1
fi

dest_dir=$1

#  Get the IDs for samples we are using in the example
spec_ids=$(cut -f 5 sample_selection/samples.manifest | cut -d "_" -f 1 | uniq)

#  Pull those samples from synapse
for specimen in $spec_ids; do
    synapse get --downloadLocation $dest_dir -q "select * from syn8408214 where specimenID='$specimen' AND fileFormat='fastq'"
done

#  Create the 'samples.manifest' file for those who want to run SPEAQeasy
Rscript pull_data/make_manifest.R -d "$dest_dir"
