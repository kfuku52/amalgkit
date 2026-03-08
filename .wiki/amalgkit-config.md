## Overview
The `amalgkit config` command helps you obtain templates for config files, which can be used as inputs for the [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select) function, together with the `metadata.tsv` file obtained from [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata). You can use the `--config` option to select a curated and pre-filled set of config files (such as `base`, `test`, `plantae`, or `vertebrate`) that you might want to tailor to your specific query. Choosing `base` creates a nearly empty set of config files, with only the formatting instructions for each file included.

## Example command
```
amalgkit config --config base
```

## Example output: `group_attribute.config`
If you would like to merge certain columns into another particular column in the `metadata.tsv` file, please specify it in this file. With this example, the `arrayexpress-developmentalstage` column is merged into the `age` column, and the `infection` column is merged into the `treatment` column. This functionality helps you browse SRA metadata that may store similar information under different attributes. Note that the column merging is done before applying `exclude_keyword.config` and `control_term.config`.

```
# case-insensitive
# regular expressions allowed
# (aggregated to)[TAB](aggregated from)

"age"	"arrayexpress-developmentalstage"
"treatment"	"infection"
```

## Example output: `control_term.config`
[`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select) uses this config file to identify "control" samples and exclude other samples from the same research, as identified by the BioProject IDs. This ensures that only non-stressed samples are included. In this example, non-control samples are removed if "mock" is found in the `treatment` column or "0 days post infection" is found in the `exp_title` column.

```
# case-insensitive
# regular expressions allowed
# (attribute)[TAB](control term)

"treatment"	"mock"
"exp_title"	"0 days post infection"
```

## Example output: `exclude_keyword.config`
This config file can be used to exclude samples containing specific bad keywords in their metadata. With this example file, any samples containing the keyword "single cell" in any of the specified columns (`exp_title`, `study_title`, `design`, `sample_title`, `sample_description`, `lib_name`, `experiment`, `treatment`, `protocol`, and `age`) are excluded, and the reason for exclusion, `single_cell`, is recorded in the `exclusion` column. Similarly, samples containing the keyword "miRNA" in either the `treatment` or `protocol` column are excluded, with the `exclusion` value marked as "small_RNA."

```
# case-insensitive
# regular expressions allowed
# (comma-delimited attributes to be searched)[TAB](arbitrary reason of exclusion)[TAB](bad keyword)

"exp_title,study_title,design,sample_title,sample_description,lib_name,experiment,treatment,protocol,age"	"single_cell"	"single cell"
"treatment,protocol"	"small_RNA"	"miRNA"
```
