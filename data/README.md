# Data provenance and local inputs

Data source: [US Accidents, Sobhan Moosavi](https://www.kaggle.com/datasets/sobhanmoosavi/us-accidents), March 2023 release. The academic analyses use New Jersey and multistate subsets of that release, not the entire national dataset.

The dataset provider lists CC BY-NC-SA 4.0 / research and non-commercial use conditions. Obtain the data from the provider and check the current terms before reuse. Dataset rights are separate from the repository's code. No dataset license is replaced by this repository.

CSV, RDS, text corpora and original course submission PDFs are deliberately excluded from the public repository. Original PDF covers contain personal identifiers. This repository keeps method source and aggregate evidence, without those covers or administrative role-assignment pages.

For team members with the original D3 and D4 ZIP archives, `tools/restore_data.py` restores provided data locally. Other readers can inspect all method source and aggregate results, but must obtain or regenerate compatible inputs before running analyses. See [reproduction notes](../docs/reproduce.md).

`Severity` measures traffic impact/disruption on a scale of 1–4. It does not measure injury or fatality. Geographic frequency and the share of records with Severity ≥ 3 must not be interpreted as exposure-adjusted road risk.
