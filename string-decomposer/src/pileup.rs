//! Get pileup for a unit from encodings.
//! # Example
//! ```rust
//! # use string_decomposer::*;
//! # use string_decomposer::pileup::*;
//! let seq = b"TTTACTACTACTGCGCGCTTT";
//! let units = vec![b"ACT".as_slice(), b"GC".as_slice()];
//! let param = StringDecompParameters::new(1,-1,-1,-1);
//! let decomposed = DecomposedSeq::new(seq, &units, &param);
//! let encodes_of_first:Vec<_> = decomposed.encoding().iter().filter(|e|e.id() == 0).collect();
//! let pileup = pileup(&encodes_of_first, &units[0]);
//! ```
use super::Encoding;

const BASES: &[u8; 4] = b"ACGT";

const fn base_table() -> [usize; 256] {
    let mut table = [0; 256];
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'a' as usize] = 0;
    table[b'c' as usize] = 1;
    table[b'g' as usize] = 2;
    table[b't' as usize] = 3;
    table
}
const BASE_TABLE: [usize; 256] = base_table();

/// Return the pileup of the encodings. [i][j] returns the fraction of the j-th base at the i-th position.
/// j: 0 <-> A,
///  : 1 <-> C,
///  : 2 <-> G,
///  : 3 <-> T.
pub fn pileup<E: std::borrow::Borrow<Encoding>>(encodings: &[E], unit: &[u8]) -> Vec<Vec<f64>> {
    if unit.iter().any(|b| !BASES.contains(b)) {
        panic!("Please remove bases other than ACGT from the unit definition.");
    }
    let id = encodings[0].borrow().id();
    if encodings.iter().any(|e| e.borrow().id() != id) {
        for enc in encodings.iter().map(|e| e.borrow()) {
            eprintln!("{enc:?}");
        }
        panic!("Different units shouldn't be accepted.");
    }
    // A,C,G,T
    let mut fractions = vec![vec![0f64; 4]; unit.len()];
    for encoding in encodings.iter().map(|e| e.borrow()) {
        let mut u_pos = 0;
        use super::Op;
        for op in encoding.ops() {
            match op {
                Op::Del(_) => u_pos += 1,
                Op::Ins(_) => {}
                &Op::Match(r_base, _u_base) => {
                    let r_base = BASE_TABLE[r_base as usize];
                    fractions[u_pos][r_base] += 1f64;
                    u_pos += 1;
                }
            }
        }
    }
    for column in fractions.iter_mut() {
        let sum: f64 = column.iter().sum();
        column.iter_mut().for_each(|x| *x /= sum);
    }
    fractions
}
