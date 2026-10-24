# Platform scope

ZRust 0.1.0 is qualified for source use on Apple Silicon macOS, Rust target `aarch64-apple-darwin`. Rust/Cargo 1.98.0 was tested. A lower minimum toolchain version has not been established. No binary distribution is provided.

Scientific calculations remain single-threaded. CPU workers prepare completed result rows. Eight formatting workers are the qualified default and the maximum used. Small outputs remain serial. Other platforms and universal speedups are not established.

The SUM149PT olaparib screen and hTERT-RPE1 olaparib screen matched across tested result fields against pinned Python DrugZ. In the HeLa olaparib screen, two genes with equal scores appeared in the opposite adjacent rank order. Statistical results and significant-gene membership were unchanged. This comparison concerns the pinned implementation, not recreation of original-publication tables.

See the [manual](ZRust-User-Manual.pdf) for operation and the [manuscript](../manuscript/ZRust-Manuscript.pdf) for methods, results and limitations. Source compilation on another platform does not establish scientific support there.
