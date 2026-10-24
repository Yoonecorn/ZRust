//! Pinned pandas 1.5.3 / NumPy 1.23.5 indirect introsort permutation.
//! Equal keys affect variance-window membership: do not substitute a host sort.
//! Algorithm provenance and BSD notice: tie-order-requalification/reconstruction.

pub(crate) fn descending(n: usize, less: impl Fn(usize, usize) -> bool) -> Vec<usize> {
    // pandas reverses values AND source indices before ascending argsort.
    let mut a: Vec<_> = (0..n).rev().collect();
    if n < 2 {
        return a;
    }
    let depth = 2 * (usize::BITS - 1 - n.leading_zeros()) as i32;
    let mut stack = vec![(0, n - 1, depth)];
    while let Some((mut lo, mut hi, mut depth)) = stack.pop() {
        // NumPy checks depth only when entering a popped range, not per split.
        if depth < 0 {
            heap(&mut a[lo..=hi], &less);
            continue;
        }
        while hi - lo > 15 {
            let mid = lo + (hi - lo) / 2;
            if less(a[mid], a[lo]) {
                a.swap(mid, lo);
            }
            if less(a[hi], a[mid]) {
                a.swap(hi, mid);
            }
            if less(a[mid], a[lo]) {
                a.swap(mid, lo);
            }
            let pivot = a[mid];
            let mut i = lo;
            let mut j = hi - 1;
            a.swap(mid, j);
            loop {
                i += 1;
                while less(a[i], pivot) {
                    i += 1;
                }
                j -= 1;
                while less(pivot, a[j]) {
                    j -= 1;
                }
                if i >= j {
                    break;
                }
                a.swap(i, j);
            }
            a.swap(i, hi - 1);
            depth -= 1;
            // Median sentinels guarantee lo<i<hi, so these bounds cannot wrap.
            if i - lo < hi - i {
                stack.push((i + 1, hi, depth));
                hi = i - 1;
            } else {
                stack.push((lo, i - 1, depth));
                lo = i + 1;
            }
        }
        for i in lo + 1..=hi {
            let saved = a[i];
            let mut j = i;
            while j > lo && less(saved, a[j - 1]) {
                a[j] = a[j - 1];
                j -= 1;
            }
            a[j] = saved;
        }
    }
    a.reverse();
    a
}

fn heap(a: &mut [usize], less: &impl Fn(usize, usize) -> bool) {
    // Local one-based heap indices mirror NumPy, with checked slice indexing.
    fn sift(
        a: &mut [usize],
        root: usize,
        n: usize,
        saved: usize,
        less: &impl Fn(usize, usize) -> bool,
    ) {
        let mut i = root;
        // Check parent bound before multiplication, avoiding index overflow.
        while i <= n / 2 {
            let mut child = 2 * i;
            if child < n && less(a[child - 1], a[child]) {
                child += 1;
            }
            if !less(saved, a[child - 1]) {
                break;
            }
            a[i - 1] = a[child - 1];
            i = child;
        }
        a[i - 1] = saved;
    }
    let mut n = a.len();
    for root in (1..=n / 2).rev() {
        let saved = a[root - 1];
        sift(a, root, n, saved, less);
    }
    while n > 1 {
        let saved = a[n - 1];
        a[n - 1] = a[0];
        n -= 1;
        sift(a, 1, n, saved, less);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_proven_order_arrays() {
        let fixtures = include_str!("../tie-order-requalification/fixtures/order-regressions.tsv");
        for line in fixtures.lines() {
            let fields: Vec<_> = line.split('\t').collect();
            let expected: Vec<usize> = fields[3]
                .split(',')
                .filter(|x| !x.is_empty())
                .map(|x| x.parse().unwrap())
                .collect();
            let actual = if fields[1].starts_with("float") {
                let v: Vec<f64> = fields[2]
                    .split(',')
                    .filter(|x| !x.is_empty())
                    .map(|x| x.parse().unwrap())
                    .collect();
                descending(v.len(), |i, j| v[i] < v[j])
            } else {
                let v: Vec<u64> = fields[2]
                    .split(',')
                    .filter(|x| !x.is_empty())
                    .map(|x| x.parse().unwrap())
                    .collect();
                descending(v.len(), |i, j| v[i] < v[j])
            };
            assert_eq!(actual, expected, "{}", fields[0]);
        }
    }
}
