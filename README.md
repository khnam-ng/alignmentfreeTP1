# RENDU - PHYG_KMERS_Phase2

## Khanh Nam NGUYEN

#### a few minutes late cause of new data folder. Its too large and I can't push it to github then I have to figure out how to reset HEAD and rebase the branch.

### Rendu_1: At this point, I did a mistake in nested loops where the execution time has been doubled.
real	8m35.929s
user	8m29.880s
sys	0m6.122s

### Rendu_1_correction: Fast and Furious
real	0m54.228s
user	0m53.564s
sys	0m0.666s

### Rendu_1_correction + heapq + sketch_size = 1000: try to implement sketch with size = 1000 and heapq by using max heap. Fast and Furious 2.
GCA_000069965.1 GCA_000008865.2 0.0
GCA_000069965.1 GCA_000005845.2 0.0
GCA_000069965.1 GCA_030271835.1 0.014541765302671627
GCA_000069965.1 GCA_000013265.1 0.00025006251562890725
GCA_000008865.2 GCA_000005845.2 0.11607142857142858
GCA_000008865.2 GCA_030271835.1 0.0
GCA_000008865.2 GCA_000013265.1 0.1001100110011001
GCA_000005845.2 GCA_030271835.1 0.0
GCA_000005845.2 GCA_000013265.1 0.1714174150722374
GCA_030271835.1 GCA_000013265.1 0.0

real	0m12.456s
user	0m12.377s
sys	0m0.061s

### Rendu_1_correction + heapq + sketch_size = 10000: increase the size of sketch to 10000. Of course this action increased also the accuracy of Jaccard similarity
Computing Jaccard similarity for all pairs of samples
GCA_000069965.1 GCA_000008865.2 0.00016002560409665546
GCA_000069965.1 GCA_000005845.2 0.00023338779048444637
GCA_000069965.1 GCA_030271835.1 0.01708706265256306
GCA_000069965.1 GCA_000013265.1 0.00020004000800160032
GCA_000008865.2 GCA_000005845.2 0.1397959765202029
GCA_000008865.2 GCA_030271835.1 0.00020004000800160032
GCA_000008865.2 GCA_000013265.1 0.11388344323649974
GCA_000005845.2 GCA_030271835.1 0.0004502025911660247
GCA_000005845.2 GCA_000013265.1 0.1916583912611718
GCA_030271835.1 GCA_000013265.1 0.00030009002700810244

real	0m12.905s
user	0m12.855s
sys	0m0.048s

### Try to implement HyperMinHash but did not success
CM000679.2 CM000678.2 1.0
CM000679.2 CM055491.2 1.0
CM000679.2 CM001012.3 0.0
CM000679.2 CM055492.2 1.0
CM000679.2 CM001004.3 1.0
CM000678.2 CM055491.2 1.0
CM000678.2 CM001012.3 0.0
CM000678.2 CM055492.2 1.0
CM000678.2 CM001004.3 1.0
CM055491.2 CM001012.3 0.0
CM055491.2 CM055492.2 1.0
CM055491.2 CM001004.3 1.0
CM001012.3 CM055492.2 0.0
CM001012.3 CM001004.3 0.0
CM055492.2 CM001004.3 1.0

real	4m39.455s
user	4m38.402s
sys	0m0.913s

