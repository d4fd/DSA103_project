Pipeline

1. Activate the virtual environment (`.venv`), install all packages listed in `requirements.txt`, and ensure **DSA103_PROJECT** is the working directory.
2. Execute the pipeline with:
   `python scripts/run.py`
3. Output files are written to `data/processed/` in `.csv` format.

Additional notes

* The original R script is still included in the repository.
* Helper functions that were not explicitly defined in the original R workflow have been re-implemented in Python and are available in `scripts/helper_functions.py`.
* `test.ipynb` can be used for testing and visualization.

Repository structure

```
.
├── README.md
├── data
│   ├── figures
│   │   ├── clustermap_molecularDescriptors.png
│   │   └── correlation_matrix_molecularDescriptors.png
│   ├── processed
│   │   ├── chemistry.csv
│   │   ├── compounds.csv
│   │   └── presAbs.csv
│   └── raw_data
│       └── mtbs_tropical_annotations.tsv
├── docs
│   ├── DSA103-Lecture21-VisualizationReview,ProjectStart.pdf
│   ├── DataStoryTemplate.docx
│   ├── DocumentationDataManagementTemplate.docx
│   ├── ReadingScientificPapers.docx
│   ├── Walker et al. - 2023 - Leaf metabolic traits reveal hidden dimensions of plant form and function.pdf
│   └── permissions
│       ├── Data_Licence
│       └── Scripts_Licence.txt
├── requirements.txt
└── scripts
    ├── derive_chemistry.R
    ├── derive_chemistry.py
    ├── helper_functions.py
    ├── run.py
    └── test.ipynb
```
