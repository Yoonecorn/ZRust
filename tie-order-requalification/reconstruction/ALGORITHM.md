# Explicit pinned local-window permutation

Reconstruction authority: installed pandas 1.5.3 `core/sorting.py:nargsort` and
NumPy v1.23.5 `aquicksort_` and `aheapsort_`; frozen source/notice is adjacent.
This specification concerns finite local-window keys, not gene-ranking ties.

The output is a permutation of retained source positions. Start the index
array with n-1, n-2, ..., 0. Indirectly sort these indices ascending by their
referenced numeric keys using the exact algorithm below, then reverse the
index array. This implements pandas' descending reverse-input/reverse-output
transformation. Index labels and original positions are not comparison
tie-breakers. Positions identify elements and initialize the permutation only.

Comparator is strict numeric `<`. Equal numeric values never compare less,
including -0.0 and +0.0. NaN and infinities are rejected by the existing finite
input domain before sorting. Paired keys retain exact nonnegative integer
tokens if every original value in the selected column is an i64/u64 integer;
otherwise use the already validated f64 values. Unpaired keys are existing
raw-control f64 means, with no formula change. Integer key preservation changes
only order comparisons; normalization, FC and output arithmetic remain f64.
Column dtype inference precedes exclusions, as in pandas loading. Token spelling
does not break numeric equality. This numeric-domain engine has no architecture
or host-sort dependency. Compatibility is defined by the pinned numeric domain,
not pandas object-string sorting outside valid legacy scientific execution.

For n<2 return the initial permutation. Otherwise maintain inclusive [lo,hi]
ranges and a signed depth budget, initially 2*floor(log2(n)). At the beginning
of processing a popped range only, if depth<0 heapsort that entire range.
Within a range, while hi-lo>15 (at least 17 elements):

1. Let mid=lo+floor((hi-lo)/2). Swap mid/lo if key(mid)<key(lo), then hi/mid
   if key(hi)<key(mid), then mid/lo if key(mid)<key(lo).
2. Save key(mid) as pivot; swap mid with hi-1. Set i=lo, j=hi-1.
3. Increment i at least once, then while key(i)<pivot. Decrement j at least
   once, then while pivot<key(j). If i>=j stop; otherwise swap i/j and repeat
   both scans. Equal keys stop scans and can be swapped across the partition.
4. Swap i with hi-1, placing the pivot. Decrement depth once. Push the larger
   partition with this depth and continue the smaller one. If i-lo<hi-i,
   push [i+1,hi], continue [lo,i-1]; otherwise push [lo,i-1], continue [i+1,hi].
   Do not add a negative-depth check inside this partition loop.

For the remaining <=16-element range, insertion-sort left to right: save each
index, shift preceding indices right only while saved key<preceding key, then
insert. Equal keys do not shift during insertion. Pop the last stacked range
and its stored depth; repeat until the stack is empty.

Heapsort uses local one-based positions. Build a max heap from floor(length/2)
down to 1. Choose the right child only when left key<right key; equal children
choose left. Sift a saved element only while saved key<chosen child key; equality
stops. Extract by saving the last element, moving root to the last position,
decreasing heap length and sifting the saved last element from root. Continue
until one element remains. Preserve these exact equality decisions.

The Python executable candidate uses no NumPy/pandas operation internally.
The Rust strategy is a generic indirect permutation engine using explicit
swaps, scans, stack and heap operations, initialized in reverse order. All
indexing is bounds checked, and partition invariants prevent underflow.
Sorting does not modify counts, guide labels or input ordinals.

Ties affect the membership of 2H variance windows and the assignment of step-25
SD blocks, including the final-tail overwrite; they can therefore change guide
Z, gene scores and FDR. Original G3 correctly failed under the former stable
tie policy. Only a successful full-order proof authorizes implementation.
