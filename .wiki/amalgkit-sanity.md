## Overview
`amalgkit sanity` scans the working directory to verify the presence of all outputs for both quant and `amalgkit getfastq` for every sample listed in `metadata.tsv`. It also checks for the requisite index files for quantification. If any discrepancies are detected, such as a missing expression file for a sample from `amalgkit quant`, `amalgkit sanity` promptly notifies the user. This function proves particularly beneficial when handling vast numbers of samples, as manual verification of all outputs becomes cumbersome.

## Usage
To check for the presence of FASTA index files, outputs from `amalgkit quant`, and outputs from `amalgkit getfastq`, use the following command. If certain checks are not required, simply omit the corresponding arguments.

```
amalgkit sanity --metadata /PATH/TO/metadata.tsv --out_dir /WORKING/DIRECTORY/ --index --quant --getfastq
```

For a comprehensive check, you can use the --all option:

```
amalgkit sanity --metadata /PATH/TO/metadata.tsv --out_dir /WORKING/DIRECTORY/ --all
```

STDOUT will display information about any missing files and suggest appropriate commands to rerun the missing samples.
