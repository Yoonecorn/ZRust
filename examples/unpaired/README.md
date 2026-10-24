# Unequal-arm unpaired teaching analysis

The input is an unchanged copy of the included 80-guide, 20-gene synthetic fixture. It is not a biological experiment. Sample columns are C1, C2, T1, T2. Run from the repository root after building. Use a fresh output folder.

<!-- command: example-unpaired -->
```sh
mkdir results-example-unpaired
target/release/zrust -i examples/unpaired/counts.tsv -c C1,C2 -x T1 -unpaired --minobs 1 \
  --half_window_size 5 -q \
  -o results-example-unpaired/genes.tsv -f results-example-unpaired/guides.tsv
cmp results-example-unpaired/genes.tsv examples/unpaired/expected/genes.tsv
cmp results-example-unpaired/guides.tsv examples/unpaired/expected/guides.tsv
```

Expected result: two complete files. Both comparisons print nothing and return success. The first gene is G18. Gene and guide columns are explained in the [user manual](../../docs/ZRust-User-Manual.pdf).

Controls C1 and C2 are compared with treatment T1 as group means. The unused T2 input column is intentional.
