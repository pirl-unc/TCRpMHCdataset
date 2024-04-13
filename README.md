[![Tests](https://github.com/pirl-unc/TCRpMHCdataset/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/TCRpMHCdataset/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/pirl-unc/TCRpMHCdataset/badge.svg?branch=main)](https://coveralls.io/github/pirl-unc/TCRpMHCdataset?branch=main)
<a href="https://pypi.python.org/pypi/tcrpmhcdataset/">
<img src="https://img.shields.io/pypi/v/tcrpmhcdataset.svg?maxAge=1000" alt="PyPI" />
</a>


# TCRpMHCdataset

TCRpMHCdataset is a library to help handle the many-to-many nature of TCR::pMHC data in a sanity preserving manner. It loads tabular data (.CSV) and converts it into a dataset object that can be indexed to get pairs of TCR and pMHC objects that retain a list of all other cognate -optes. The dataset is designed with the flexibility to be used with minimal changes for most TCR:pMHC related tasks. It comes packaged with a notion of directionality: TCR -> pMHC (De-Orphanization) or pMHC -> TCR (TCR design). For training and evaluating machine learning models, the dataset object implements a `split` function which robustly handles the splitting of the data into balanced train and test splits that can be stratified by the allele frequences. These splits can also explicitly hold out epitopes, epitope::allele combinations, and even TCRs that are present in the test set from the training set. 

## Installation

----------------------------------------------------------------------

### From PyPI
```bash
pip install tcrpmhcdataset
```

### From source

Clone this repository into your working directory and run the standard installation command:
```bash
git clone https://github.com/pirl-unc/TCRpMHCdataset.git
cd TCRpMHCdataset
pip install . 
```

## Usage Examples

----------------------------------------------------------------------

### Loading a Dataset
```python
from tcrpmhcdataset import TCRpMHCdataset
# Define a TCR -> pMHC dataset
deorph_dataset = TCRpMHCdataset(source='tcr', target='pmhc', use_mhc=False, use_pseudo=True, use_cdr3=True, use_both_chains=False)
deorph_dataset.load('test_data/sampled_paired_data_cleaned.csv')

In [1]: print(deorph_dataset)
Out[1]: 'TCR:pMHC Dataset of N=6833. Mode:tcr -> pmhc.'
```

### Indexing an Item
```python
from tcrpmhcdataset import TCRpMHCdataset
# Define a pMHC -> tcr dataset
design_dataset = TCRpMHCdataset(source='pmhc', target='tcr', use_mhc=False, use_pseudo=True, use_cdr3=True, use_both_chains=False)
design_dataset.load('test_data/sampled_paired_data_cleaned.csv')
pmhc, tcr = design_dataset[0]

In [1]: pmhc
Out[1]: pMHC(peptide="LIDFYLCFL", hla_allele="HLA-A*02:01", reference={'MIRA:eEE226', 'MIRA:eEE240', 'MIRA:eOX54', 'MIRA:eEE224', 'MIRA:eXL37', 'MIRA:eOX52', 'MIRA:eOX43', 'MIRA:ePD76', 'MIRA:eHO130', 'MIRA:eQD137', 'MIRA:eXL31', 'MIRA:eHH175', 'MIRA:eOX56', 'MIRA:eXL30', 'MIRA:eXL27'}, use_pseudo=True, use_mhc=False)

In [2]: tcr
Out[2]: TCR(cdr3a="None", cdr3b="CSAQDRTSNEQFF",
                trav="None", trbv="TRBV20-1", 
                traj="None", trbj="TRBJ2-1",
                trad="None", trbd="None",
                tcra_full="None", tcrb_full="MLLLLLLLGPGISLLLPGSLAGSGLGAVVSQHPSWVICKSGTSVKIECRSLDFQATTMFWYRQFPKQSLMLMATSNEGSKATYEQGVEKDKFLINHASLTLSTLTVTSAHPEDSSFYICSAQDRTSNEQFFGPGTRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG",
                reference={'MIRA:eEE226'}, use_cdr3b=True)
```

### Splitting a Dataset
```python
from tcrpmhcdataset import TCRpMHCdataset
# Define a TCR -> pMHC dataset
design_dataset = TCRpMHCdataset(source='pmhc', target='tcr', use_mhc=False, use_pseudo=True, use_cdr3=True, use_both_chains=False)
design_dataset.load('test_data/sampled_paired_data_cleaned.csv')

# Split on Epitope and then pull out a validation set
train_dataset, test_dataset = design_dataset.split(test_size=0.2, balance_on_allele=True, split_on=['Epitope'])
train_dataset, val_dataset = train_dataset.split(test_size=0.1, balance_on_allele=True, split_on=['Epitope', 'Allele'])

# Split on pMHC 
train_dataset2, test_dataset2 = design_dataset.split(test_size=0.2, balance_on_allele=True, split_on=['Epitope', 'Allele'])

# Split on TCR
train_dataset3, test_dataset3 = design_dataset.split(test_size=0.2, balance_on_allele=True, split_on=['CDR3b', 'CDR3a'])
```

## Motivation

----------------------------------------------------------------------

T-cell Receptors (TCRs) are highly specific pattern recognition receptors that allow T-cells to recognize non-self molecular motifs. Though necessarily highly specific in order to avoid self-reactivity, the phenomenon of cross-reactivity is required to maintain physiologically manageable numbers of T-cell clones while still ensuring sufficient protection (See [Why must T cells be cross-reactive](https://www.nature.com/articles/nri3279)). While fascinating from a biological perspective, handling the data for this complex many-to-many mapping is another story. This project started off as three abstractions that I wrote for [this project](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Alnulv8AAAAJ&citation_for_view=Alnulv8AAAAJ:2osOgNQ5qMEC) to help think about TCRs and pMHCs are distinct entities with attributes pertaining to what is known about them and a dataset class that could quickly load data from a tabular format, split them in a balanced and meaningful manner, and be easily indexed for model training or evaluation. After having used these abstractions in a few projects, I reckon that they could be useful to others as well and decided to package them up and release them as a library.


## Modules

----------------------------------------------------------------------

See the [documentation](https://pirl-unc.github.io/TCRpMHCdataset/tcrpmhcdataset.html) for more information.

#### TCRpMHCdataset

The `TCRpMHCdataset` object is the main object that is used to load, index, and split the data. 

#### TCR

The TCR object is a frozen dataclass object that stores information such as the 'CDR3b' sequence, 'CDR3a', 'Vb', 'Jb', etc. but also includes a set of pMHCs that this TCR is reactive against as well as a set of references that support the TCR. 

#### pMHC

Like the TCR's implementation, each pMHC object similarly contains key information including the 'peptide', 'allele' but also includes a set of TCRs that are reactive against it as well as a set of references that support the pMHC. It also computes the full HLA-sequence as well as the pseudosequence and caches these for future reference. 

*MHC Sequence*

Major Histocompatibility Complex (MHC) sequences are mapped using the parsed allele level information using the [IMGT HLA]() database. MHC proteins are sometimes annotated with mutations in relation to a known allele:

- "HLA-B\*08:01 N80I mutant"

If picked up by `mhcgnomes` These are passed in to the pMHC object which makes the necessary changes to the sequence and caches the new sequence for future reference.

*Pseudo Sequence*

Pseudosequences are derived from the full MHC sequence that are predicted to be in contact with the peptide given proximity and a polymorphism based estimator (Introduced in [netMHCpan](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000796#s4)). All pseudo-sequenes are 34 Amino Acids long and are used in different immuno-informatics pipelines as a reduced representation of the full MHC seqeunce. In this package, the pseudo-sequences are mapped from allele level information, similar to MHC sequences. They are not updated given a mutation, but instead use the canonical allele's pseudo-sequence, if available.

*Allele Imputation strategy*

Given the sparsity of the data, every single datapoint is of critical importance. In the event that only the serotype information was provided (e.g. HLA-A2), an "eager" imputation strategy is included to try and impute a common allele (e.g. HLA-A\*02:01). This is done in a rudimentary manner by guessing the allele field from :01 -> :10 and seeing if there exists an MHC sequence from IMGT that matches the allele. This strategy should ONLY be used if using pseudosequence level information as the pseudosequence often is highly conserved within serotype. 

## Contributing

----------------------------------------------------------------------

Community help is always appreciated. No contribution is too small, and we especially appreciate efficiency and usability improvements such as better documentation, tutorials, tests, or code cleanup. If you're looking for a place to start, check out the issues labeled "good first issue" in the issue tracker.

#### Project scope
The `TCRpMHCdataset`, `TCR`, and `pMHC` classes are designed with object oriented principles in mind to think of these protein complexes as individual units with intra-and interdependent interactions, not rows in a dataframe. To this end, any and all efforts to expand the expressivity of these objects with available data is most certainly within scope. It is our hope to soon incorporate structure level information as well. 

All committed code to `TCRpMHCdataset` should be suitable for regular research use by practioners.

If you are contemplating a large contribution, such as the addition of a new class, modality, or data-structure altogether, please first open an issue on GH (or email us at dkarthikeyan1@unc.edu) to discuss and coordinate the work.

#### Making a contribution
All contributions can be made as pull requests on Github. One of the core developers will review your contribution. As needed the core contributors will also make releases and submit to PyPI.

A few other guidelines:

 * `TCRpMHCdataset` is written for Python3 on Linux and OS X. We can't guarantee support for Windows.
 * All functional modifications should be documented using [numpy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html) with corresponding with unit tests.
 * Please use informative commit messages.
 * Bugfixes should be accompanied with test that illustrates the bug when feasible.
 * Contributions are licensed under Apache 2.0
 * All interactions must adhere to the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/1/4/code-of-conduct/).


## References

- [Conditional Generation of Antigen Specific T-cell Receptor Sequences](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=Alnulv8AAAAJ&citation_for_view=Alnulv8AAAAJ:2osOgNQ5qMEC)
- [Why must T cells be cross-reactive](https://www.nature.com/articles/nri3279)
- [Mhcgnomes](https://github.com/pirl-unc/mhcgnomes)
- [Tidytcells](https://tidytcells.readthedocs.io/en/latest/#)
- [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/1/4/code-of-conduct/)