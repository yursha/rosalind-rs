An idiomatic, type-safe bioinformatics library and CLI toolkit written in Rust, built and validated using Rosalind challenges.

This project is under active development. The core library architecture is being built incrementally by implementing robust, production-grade solutions to the Rosalind bioinformatics problem set. 
It is not yet published on crates.io, but the library engine and existing CLI tools are fully functional and test-verified locally.

* **`src/dna`**: Nucleotide operations, transcription, reverse complements.
* **`src/rna`**: Translation, codon mapping, splice site analysis (Planned).
* **`src/protein`**: Peptide profiling, molecular weight calculations (Planned).
* **`src/bin/`**: A suite of lightweight, high-performance CLI utilities built on top of the core library.
