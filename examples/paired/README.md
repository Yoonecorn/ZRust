# Paired teaching analysis

The input is an unchanged copy of the included 80-guide, 20-gene synthetic fixture. It is not a biological experiment. Sample columns are C1, C2, T1, T2. Run from the repository root after building. Use a fresh output folder.

<!-- command: example-paired -->
```sh
mkdir results-example-paired
target/release/zrust -i examples/paired/counts.tsv -c C1,C2 -x T1,T2 \
  --half_window_size 5 -q \
  -o results-example-paired/genes.tsv -f results-example-paired/guides.tsv
cmp results-example-paired/genes.tsv examples/paired/expected/genes.tsv
cmp results-example-paired/guides.tsv examples/paired/expected/guides.tsv
```

Expected result: two complete files. Both comparisons print nothing and return success. The first gene is G18. Gene and guide columns are explained in the [user manual](../../docs/ZRust-User-Manual.pdf).

C1 pairs with T1; C2 pairs with T2.
