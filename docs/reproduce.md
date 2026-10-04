# Reproduction and execution notes

This is a curated academic research archive. Source is available for review, and several analyses can start from the submitted feature dataset. It is not a fully locked, one-command pipeline. The original experiment environment and every intermediate file were not included in the submissions.

## What was verified for this publication

- All 19 R files parse with R, without loading analysis packages or running models.
- Notebook structure and Python cell syntax were checked after removing embedded outputs and local paths.
- D4 final cluster counts were recomputed from the submitted assignment CSV: 5,195 / 4,149 / 2,476 / 2,562 / 307, totaling 14,689.
- The New Jersey feature CSV contains 8,923 records and 19 columns.
- Local ZIP restoration was exercised against both supplied archives.

These checks verify the archive and its presentation. They do not establish that every analysis reruns end to end on a fresh computer. No models were refitted for this publication.

## Data setup

Use Python 3 and R. For team members with the original submissions, run from the repository root:

```sh
python3 tools/restore_data.py --d3 /path/to/ID-3-USAccidents-D3.zip --d4 /path/to/ID-3-USAccidents-D4.zip
```

The helper copies CSV/RDS files into ignored `data/` and `workspace/` directories. It does not execute R or install packages. The public repository does not supply those course archives. Other readers should obtain the March 2023 US Accidents release from the original provider and adapt/regenerate the subsets using the published preprocessing source. Random sampling, package changes and unavailable historical intermediates can prevent exact reproduction of the archived numbers.

Each R file declares its libraries near the top. Install them in an isolated R environment before execution. Several original scripts attempt package installation automatically if packages are missing; inspect that behavior before running. Dependencies have not been pinned, and source may require adaptation to current package APIs. The Python notebook imports pandas, NumPy, matplotlib, seaborn, Plotly, kaleido, folium, branca, scikit-learn, scipy, missingno, geopandas, shapely and h3. Map tiles require internet access.

## Practical entry points

Run data-dependent R scripts from the corresponding `workspace/d3` or `workspace/d4` directory, so the published relative input paths resolve. For example:

```sh
cd workspace/d4
Rscript ../../analysis/d4/scripts/FAMD_NJ_PMAAD_variant2_noCounty_noLatLng.R
Rscript ../../analysis/d4/scripts/LDA.R
Rscript ../../analysis/d4/scripts/CA.R
Rscript ../../analysis/d4/scripts/Geo.R
Rscript "../../analysis/d4/scripts/Scripts Clustering/profiling_final_k5.R"
```

These inputs exist in the supplied D4 archive after restoration; successful runtime also requires their respective packages and compatible APIs. The profiling script reads `Datasets/10_base_final_para_profiling_k5.csv`. It automatically creates local plot and table folders.

The MCA entry point searches for the repository root and reads `data/US_Accidents_Final_FE.csv`; launch it from the repository root:

```sh
Rscript analysis/d4/scripts/acm_us_accidents_final.R
```

Open `analysis/d4/geospatial_statistics.ipynb` with a Python 3 kernel. Its input path supports launch from the notebook directory or repository root; run cells in order. Outputs are cleared in Git so notebook views do not embed large maps or environment-specific data.

D3 county temporal clustering can start from the restored feature CSV:

```sh
cd workspace/d3
Rscript "../../analysis/d3/scripts/ts_clustering.R"
```

## Pipeline branches and missing intermediates

The multistate branch starts at `preprocesado_clustering_US_accidents.R`, requiring `US_Accidents_March23.csv` in its working directory. It generates `Datasets/clustering/03_sample_15000_clean.csv` and `04_base_clustering.csv`, needed by D3 k-prototypes/CURE-like comparison and D4 multistate clustering. These processed intermediate files are absent from the submitted ZIPs; the raw 15,000-row sample is included, but it is not a drop-in replacement for the original full-dataset input.

The New Jersey annex branch follows EDA → imputation → outlier treatment → feature engineering. Individual scripts read inputs from `Datasets/` but often write to the working directory. Move generated intermediates into `Datasets/` and match filename case before running the next stage. The original combined script is preserved as an alternative historical entry point; inspect its assumptions rather than mixing both paths.

The two D4 clustering scripts require `Datasets/clustering/04_base_clustering.csv`. The exploratory script considers dendrogram heights; the final script explicitly selects five clusters. The submitted final assignment CSV supports profiling without refitting that branch. Restoring data does not reconstruct missing upstream intermediates or guarantee the same fitted partition.

## Curation changes

Only portability and presentation edits were made: private absolute paths became relative paths; placeholder input/output paths were replaced; two inconsistent RDS filename cases were normalized; profiling uses the working directory; MCA detects the public repository layout; notebook outputs and environment metadata were cleared. Notebook labels describing Severity ≥ 3 as medical severity were corrected to traffic impact. Statistical methods and archived result counts were not changed. `evidence/source-manifest.json` records original and published R source hashes.

Original submission PDFs are not published because their covers include personal identifiers and their administrative pages do not belong in the project presentation. The README and aggregate evidence provide the public technical summary.
