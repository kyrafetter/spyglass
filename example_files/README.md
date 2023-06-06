We have provided the following mini test files in order to test the installation of `spyglass`:

#### `test_peaks.bed`
A BED file of foreground genomic peak regions. These will be tested for known motif enrichment.

#### `test_background.bed`
A BED file of background genomic peak regions. These will be used as custom background, should the user choose to provide this file.

#### `test_motifs.pwm`
A PWM file to be used by `spyglass` containing three motifs to be tested for enrichment in `test_peaks.bed`.

#### `test_homer.motifs`
A motifs file to be used by `Homer` containing three motifs to be tested for enrichment in `test_peaks.bed`.

#### `test_ref.fa`
A reference sequence in FASTA format. 
