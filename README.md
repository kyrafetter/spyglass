# spyglass
spyglass is a PWM (Position Weight Matrix)-based implementaion of a DNA binding motif finder. It is a tool that identifies enriched motifs in 1D genomic regions given a list of known motifs. 

[Prerequisites](#prerequisites) | [Installation](#install) | [Basic Usage](#usage) | [spyglass Options](#options) | [File formats](#formats) | [Testing](#testing)

<a name="prerequisites"></a>
## Prerequisites

<a name="install"></a>
## Installation

<a name="usage"></a>
## Basic Usage 
The basic usage of `spyglass` is:
```
spyglass [other options] ref.fa peaks.bed motifs.pwm
```
To run `spyglass` on our mini test files (see `example-files`):
```
spyglass test_ref.fa test_peaks.bed test_motifs.pwm
```

<a name="options"></a>
## spyglass Options
`spyglass` has the following required arguments:
- `ref.fa`: reference genome in Fasta format
- `peaks.bed`: genomic regions file, for example, peak calls from ChIP-seq datasets.
-  `motifs.pwm`: PWMs of all motifs-of-interest. `spyglass` will determine whether these motifs are significantly enriched in `regions.bed`

Additionally, users may choose to specify the optional options below:
 - `-b BACKGROUND`, `--background BACKGROUND`: Use custom frequencies for A, T, C, and G. By default, all bases have a background frequency of 0.25.
 - `-o FILE`, `--output FILE`: Write output to this file. By default, output is written to stdout.

<a name="formats"></a>
## File Formats
`peaks.bed` is a tab-delimited file in BED format with no header. It contains six columns as follows:
```
chromosome start_coordinate  end_coordinate peak_ID . strand(+/-)
```
`ref.fa` is a reference genome sequence in FASTA format. It contains the name of the chromosome followed by the sequence:
```
>chr[name]
[chromosome sequence]
```

<a name="testing"></a>
## Testing

<a name="contributors"></a>
## Contributors 
