# ont simsig

Project to simulate simple Oxford Nanopore Technologies RNA signals

## signal simulation

The signal gets simulated by providing a reference sequence or generating a random reference sequence.

### segment signal simulation

We use Gaussian signal distributions provided by [ONT](data\template_median69pA.model) to simulate the signal for each 5-mer in the reference sequence.

### segment length simulation

We analyzed the segment distribution on IVT RNA data sequenced with ONT. The segment lengths were calculated by ONTs software ```nanopolish eventalign```. We fitted [negative binomial distributions](data\kmer_nbin.csv) on each 5-mer length distribution and use them to simulate segment lengths with ```ont_simig```.