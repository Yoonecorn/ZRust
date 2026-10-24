# Licensing and attribution for this source release

ZRust is an independent Rust reimplementation of DrugZ, not endorsed by the
upstream DrugZ authors. ZRust-authored material is licensed by LICENSE.
Upstream DrugZ: https://github.com/hart-lab/drugz ; reviewed commit
`eb15d34e4dd172965e618d5bb662c053066da799`.
No upstream DrugZ implementation source is included in this archive. Small
accepted outputs are retained as regression references. The upstream MIT notice
is preserved verbatim below for attribution.

## DrugZ MIT notice

```text
MIT License

Copyright (c) 2018 hart-lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Ordering algorithm provenance

The Rust permutation engine was independently reconstructed from NumPy 1.23.5
indirect argsort and pandas 1.5.3 nargsort behavior. Preserve their complete
notices, included unmodified with the algorithm specification:

- `tie-order-requalification/reconstruction/numpy-source/LICENSE.txt`.
- `tie-order-requalification/reconstruction/pandas-LICENSE.txt`.
- `tie-order-requalification/reconstruction/ALGORITHM.md`.

The NumPy notice includes its bundled-component notices in full. No upstream
source excerpts, NumPy installation or pandas installation are shipped here.

## Cargo dependencies

The direct dependencies are libm 0.2.16 (MIT) and tempfile 3.27.0
(MIT OR Apache-2.0), with default features disabled. Cargo.lock records the
exact transitive graph and registry checksums. Dependency source and binaries
are not bundled; Cargo obtains them from its configured registry/cache with
their own license notices. No vendoring or binary distribution is performed.
