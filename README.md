# Trait-based diversity with kernel hypervolumes

This repository demonstrates a trait-based diversity workflow using kernel hypervolumes in R.

The example compares the occupied trait space of two penguin species using three continuous morphological traits:
- bill length
- bill depth
- flipper length

Hypervolumes are estimated with the `hypervolume` package using Gaussian kernel density estimation. Pairwise overlap is then quantified with standard set-based metrics.

## Aim

The objective is to show a minimal but complete workflow for trait-space analysis:

1. assemble a trait matrix
2. standardise traits
3. estimate species-level hypervolumes
4. quantify overlap and uniqueness
5. visualise the resulting trait spaces

This is a methodological demonstration, not a full ecological analysis.

## Data

The example uses the `palmerpenguins` dataset, treating species as groups occupying a multidimensional trait space.

## Methods

Hypervolumes are constructed with `hypervolume_gaussian()`, which estimates a Gaussian kernel density in multidimensional space.  
Bandwidths are estimated with `estimate_bandwidth()`.  
Pairwise set operations are computed with `hypervolume_set()`, and overlap statistics are extracted with `hypervolume_overlap_statistics()`.

Traits are z-standardised before hypervolume estimation so that each axis contributes on a comparable scale.

## Requirements

Install the required packages in R:

```r
install.packages(c("hypervolume", "palmerpenguins", "dplyr", "ggplot2"))
