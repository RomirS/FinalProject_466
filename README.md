# Performance Analysis of Phylogenetic Tree Construction for SARS-CoV-2 Variants

### Running the project
1. Open `main.ipynb`
2. Start runtime in kernel
2. Run all cells

### To test a different input size
1. Set `MAX_COUNT` to be `input_size/5` because this will be the desired sequence count per variant (5 variants).
2. Set `file_loc` to `'./data/dist[X].txt'` where `[X]` is the desired input size
3. Add the desired input size to `input_lengths`
4. Rerun all cells

### Dependencies (will be installed when running notebook)
numpy, pandas, matplotlib, Bio