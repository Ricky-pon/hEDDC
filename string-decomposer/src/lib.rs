//! String decomposer -- clone
//! # Example
//! ```rust
//! # use string_decomposer::*;
//! let seq = b"TTTACTACTACTGCGCGCTTT";
//! let units = vec![b"ACT".as_slice(), b"GC".as_slice()];
//! let param = StringDecompParameters::new(1,-1,-1,-1);
//! let decomposed = DecomposedSeq::new(seq, &units, &param);
//! ```
//! # Reference
//! - original paper
//! - Morishita-sensei's implementation
//!
//!
pub mod pileup;
#[derive(Debug, Clone)]
pub struct DecomposedSeq {
    begin: usize,
    end: usize,
    score: i64,
    encodings: Vec<Encoding>,
}

#[derive(Debug, Clone)]
pub struct Encoding {
    unit_id: usize,
    ops: Vec<Op>,
}

impl Encoding {
    /// Return the index of the unit.
    pub fn id(&self) -> usize {
        self.unit_id
    }
    /// Return the operations.
    pub fn ops(&self) -> &[Op] {
        &self.ops
    }
    /// Push the unit sequence into the given buffer.
    pub fn push_read_seq(&self, buffer: &mut Vec<u8>) {
        for op in self.ops.iter() {
            match *op {
                Op::Del(_) => {}
                Op::Ins(base) => buffer.push(base),
                Op::Match(read, _) => buffer.push(read),
            }
        }
    }
}

#[derive(Clone)]
pub enum Op {
    /// Insertion to the unit.
    Ins(u8),
    /// Deletion from the unit.
    Del(u8),
    /// Mismatch between the read and the unit. `(read_base, unit_base)`
    Match(u8, u8),
}

impl std::fmt::Debug for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::Ins(arg0) => write!(f, "I({})", arg0 as char),
            Self::Del(arg0) => write!(f, "D({})", arg0 as char),
            Self::Match(arg0, arg1) if arg0 != arg1 => {
                write!(f, "X({}{})", arg0 as char, arg1 as char)
            }
            Self::Match(arg0, arg1) => write!(f, "=({}{})", arg0 as char, arg1 as char),
        }
    }
}

pub struct StringDecompParameters {
    match_score: i64,
    mismatch_score: i64,
    // Insertion to the *unit*. Should be negative.
    ins_score: i64,
    // Deletion from the *unit*. Should be negative.
    del_score: i64,
}

impl StringDecompParameters {
    pub const fn new(mat: i64, mism: i64, ins: i64, del: i64) -> Self {
        Self {
            match_score: mat,
            mismatch_score: mism,
            ins_score: ins,
            del_score: del,
        }
    }
}
impl DecomposedSeq {
    /// Return the 0-based index of the alignment region.
    /// For example, if the returned value if `(s,e)`, `&read[s..e]` is the region to be aligned.
    pub fn location(&self) -> (usize, usize) {
        (self.begin, self.end)
    }
    /// Return the encoding by the String decomposer algorithm.
    pub fn encoding(&self) -> &[Encoding] {
        &self.encodings
    }
    pub fn score(&self) -> i64 {
        self.score
    }
    /// The string decomposer algrithm.
    pub fn new<U: std::borrow::Borrow<[u8]>>(
        seq: &[u8],
        units: &[U],
        param: &StringDecompParameters,
    ) -> Self {
        let (score, path) = string_decomposer(seq, units, param);
        let begin = path.first().map(|x| x.0).unwrap_or(0);
        let end = path.last().map(|x| x.0).unwrap_or(seq.len());
        let encodings = split_into_encoding(seq, units, &path);
        Self {
            begin,
            end,
            score,
            encodings,
        }
    }
}

fn split_into_encoding<U: std::borrow::Borrow<[u8]>>(
    seq: &[u8],
    units: &[U],
    path: &[(usize, usize, usize)],
) -> Vec<Encoding> {
    let units: Vec<_> = units.iter().map(|u| u.borrow()).collect();
    let mut encodings = vec![];
    let mut buffer = vec![];
    let (mut r_pos, mut u_id, mut u_pos) = match path.first() {
        Some(res) => *res,
        None => return vec![],
    };
    for &(c_r_pos, c_u_id, c_u_pos) in path.iter().skip(1) {
        // Check if we entered the next iteration...
        if c_u_pos == 0 && u_pos != 0 {
            let enc = Encoding {
                unit_id: u_id,
                ops: buffer.clone(),
            };
            encodings.push(enc);
            buffer.clear();
            u_pos = 0;
            continue;
        }
        // Determine the operation.
        let op = match (c_r_pos - r_pos, c_u_pos - u_pos) {
            (1, 1) => Op::Match(seq[c_r_pos - 1], units[c_u_id][c_u_pos - 1]),
            (0, 1) => Op::Del(units[c_u_id][c_u_pos - 1]),
            (1, 0) => Op::Ins(seq[c_r_pos - 1]),
            _ => panic!("{},{},{},{}", c_r_pos, r_pos, c_u_pos, u_pos),
        };
        buffer.push(op);
        u_id = c_u_id;
        r_pos = c_r_pos;
        u_pos = c_u_pos;
    }
    if !buffer.is_empty() {
        let enc = Encoding {
            unit_id: u_id,
            ops: buffer,
        };
        encodings.push(enc);
    }
    encodings
}

// The first element, each "row" of the DP-table
type Row = (i64, Vec<Vec<i64>>);
#[derive(Debug, Clone)]
struct StrDecompDP {
    inner: Vec<Row>,
}

impl StrDecompDP {
    fn new(seq: &[u8], units: &[&[u8]]) -> Self {
        let row: Vec<_> = units.iter().map(|u| vec![0; u.len()]).collect();
        let inner = vec![(0, row.clone()); seq.len() + 1];
        Self { inner }
    }
    fn get_row_mut(&mut self, i: usize) -> (&mut i64, &mut [Vec<i64>]) {
        self.inner
            .get_mut(i)
            .map(|(x, y)| (x, y.as_mut_slice()))
            .unwrap()
    }
    fn get_mut(&mut self, read_pos: usize, unit_id: usize, unit_pos: usize) -> &mut i64 {
        match unit_pos {
            0 => &mut self.inner[read_pos].0,
            _ => &mut (self.inner[read_pos].1)[unit_id][unit_pos - 1],
        }
    }
    fn get(&self, read_pos: usize, unit_id: usize, unit_pos: usize) -> i64 {
        match unit_pos {
            0 => self.inner[read_pos].0,
            _ => self.inner[read_pos].1[unit_id][unit_pos - 1],
        }
    }
    fn last_max_of(&self, r_pos: usize) -> (i64, usize, usize) {
        self.inner[r_pos]
            .1
            .iter()
            .enumerate()
            .map(|(u_id, row)| (u_id, row.iter().enumerate().last()))
            .filter_map(|(u_id, x)| x.map(|(u_pos, score)| (*score, u_id, u_pos + 1)))
            .max_by_key(|x| x.0)
            .unwrap()
    }
}

// Do global  alignment...
fn string_decomposer<U: std::borrow::Borrow<[u8]>>(
    seq: &[u8],
    units: &[U],
    param: &StringDecompParameters,
) -> (i64, Vec<(usize, usize, usize)>) {
    let units: Vec<_> = units.iter().map(|u| u.borrow()).collect();
    // Initialization.
    let mut dp_table = StrDecompDP::new(seq, &units);
    {
        let (_, rows) = dp_table.get_row_mut(0);
        for row in rows.iter_mut() {
            for (i, slot) in row.iter_mut().enumerate() {
                *slot = param.del_score * (i as i64 + 1);
            }
        }
    }
    // DP
    for (i, s_base) in seq.iter().enumerate().map(|(i, b)| (i + 1, b)) {
        for (u_id, unit) in units.iter().enumerate() {
            for (j, u_base) in unit.iter().enumerate().map(|(j, b)| (j + 1, b)) {
                let mat_score = match s_base == u_base {
                    true => param.match_score,
                    false => param.mismatch_score,
                };
                let mat_score = dp_table.get(i - 1, u_id, j - 1) + mat_score;
                let del_score = dp_table.get(i, u_id, j - 1) + param.del_score;
                let ins_score = dp_table.get(i - 1, u_id, j) + param.ins_score;
                *dp_table.get_mut(i, u_id, j) = mat_score.max(del_score).max(ins_score);
            }
        }
        let (wrap, _u_id, _u_pos) = dp_table.last_max_of(i);
        *dp_table.get_mut(i, 0, 0) = wrap.max(0);
    }
    // Traceback
    let mut alnpath = vec![];
    let (mut score, mut r_pos, mut u_id, mut u_pos) = (0..seq.len() + 1)
        .map(|r_pos| {
            let (score, u_id, u_pos) = dp_table.last_max_of(r_pos);
            (score, r_pos, u_id, u_pos)
        })
        .max_by_key(|x| x.0)
        .unwrap();
    // eprintln!("{score}");
    let opt_score = score;
    while !(score == 0 && u_pos == 0) && (0 < r_pos) {
        // eprintln!("{r_pos},{u_id},{u_pos}");
        alnpath.push((r_pos, u_id, u_pos));
        let u_base = units[u_id][u_pos - 1];
        let r_base = seq[r_pos - 1];
        let mat_score = if u_base == r_base {
            param.match_score
        } else {
            param.mismatch_score
        };
        let mat_score = dp_table.get(r_pos - 1, u_id, u_pos - 1) + mat_score;
        let del_score = dp_table.get(r_pos, u_id, u_pos - 1) + param.del_score;
        let ins_score = dp_table.get(r_pos - 1, u_id, u_pos) + param.ins_score;
        if score == mat_score {
            r_pos -= 1;
            u_pos -= 1;
        } else if score == del_score {
            u_pos -= 1;
        } else if score == ins_score {
            r_pos -= 1;
        } else {
            unreachable!()
        }
        score = dp_table.get(r_pos, u_id, u_pos);
        // If reached the start of the unit, wrap around.
        if u_pos == 0 && r_pos != 0 && score != 0 {
            alnpath.push((r_pos, u_id, u_pos));
            let (prev_max, next_u_id, next_u_pos) = dp_table.last_max_of(r_pos);
            assert_eq!(prev_max, score);
            u_id = next_u_id;
            u_pos = next_u_pos;
        }
    }
    while 0 < u_pos {
        alnpath.push((r_pos, u_id, u_pos));
        u_pos -= 1;
    }
    alnpath.push((r_pos, u_id, u_pos));
    alnpath.reverse();
    (opt_score, alnpath)
}

#[cfg(test)]
mod tests {
    use super::*;
    const PARAM: StringDecompParameters = StringDecompParameters::new(1, -1, -1, -1);
    #[test]
    fn string_decomposer_test_7() {
        let read = b"TTCATTTC";
        let units = vec![b"ATTTC".as_slice()];
        let (score, ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, 6);
        let path = vec![
            (0, 0, 0),
            (0, 0, 1),
            (0, 0, 2),
            (1, 0, 3),
            (2, 0, 4),
            (3, 0, 5),
            (3, 0, 0),
            (4, 0, 1),
            (5, 0, 2),
            (6, 0, 3),
            (7, 0, 4),
            (8, 0, 5),
        ];
        assert_eq!(ops, path);
    }
    #[test]
    fn string_decomposer_test() {
        let read = b"ACGACGACG";
        let units = vec![b"ACG".as_slice()];
        let (score, ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, 9);
        let path = vec![
            (0, 0, 0),
            (1, 0, 1),
            (2, 0, 2),
            (3, 0, 3),
            (3, 0, 0),
            (4, 0, 1),
            (5, 0, 2),
            (6, 0, 3),
            (6, 0, 0),
            (7, 0, 1),
            (8, 0, 2),
            (9, 0, 3),
        ];
        assert_eq!(ops, path);
    }
    #[test]
    fn string_decomposer_test_2() {
        let read = b"TTTACGCCGACGTTT";
        let units = vec![b"ACG".as_slice(), b"CCG".as_slice()];
        let (score, ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, 9);
        let path = vec![
            (3, 0, 0),
            (4, 0, 1),
            (5, 0, 2),
            (6, 0, 3),
            (6, 1, 0),
            (7, 1, 1),
            (8, 1, 2),
            (9, 1, 3),
            (9, 0, 0),
            (10, 0, 1),
            (11, 0, 2),
            (12, 0, 3),
        ];
        assert_eq!(path, ops);
    }
    #[test]
    fn string_decomposer_test_3() {
        let read = b"TTTACGCCGACGTTT";
        let units = vec![
            b"A".as_slice(),
            b"C".as_slice(),
            b"G".as_slice(),
            b"T".as_slice(),
        ];
        let (score, _ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, read.len() as i64);
    }
    #[test]
    fn string_decomposer_test_4() {
        let read = b"TTTACGCTGACGTTT";
        let units = vec![b"ACG".as_slice(), b"CCG".as_slice()];
        let (score, ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, 7);
        let path = vec![
            (3, 0, 0),
            (4, 0, 1),
            (5, 0, 2),
            (6, 0, 3),
            (6, 1, 0),
            (7, 1, 1),
            (8, 1, 2),
            (9, 1, 3),
            (9, 0, 0),
            (10, 0, 1),
            (11, 0, 2),
            (12, 0, 3),
        ];
        assert_eq!(path, ops);
        let read = b"TTTACGCCTGACGTTT";
        let units = vec![b"ACG".as_slice(), b"CCG".as_slice()];
        let (score, ops) = string_decomposer(read, &units, &PARAM);
        assert_eq!(score, 8);
        let path = vec![
            (3, 0, 0),
            (4, 0, 1),
            (5, 0, 2),
            (6, 0, 3),
            (6, 1, 0),
            (7, 1, 1),
            (8, 1, 2),
            (9, 1, 2),
            (10, 1, 3),
            (10, 0, 0),
            (11, 0, 1),
            (12, 0, 2),
            (13, 0, 3),
        ];
        assert_eq!(path, ops);
    }
    #[test]
    fn string_decomposer_test_5() {
        let read = b"AATAA";
        let units = vec![b"A".as_slice(), b"CCG".as_slice()];
        let decomp = DecomposedSeq::new(read, &units, &PARAM);
        assert_eq!(decomp.score, 3);
        eprintln!("{:?}", decomp.encoding());
        assert_eq!(decomp.encoding().len(), 5);
    }
    #[test]
    fn string_decomposer_test_6() {
        let read = b"AACGTCGTCAAGTAAA";
        let units = vec![b"CGT".as_slice()];
        let decomp = DecomposedSeq::new(read, &units, &PARAM);
        assert_eq!(decomp.score, 7);
        assert_eq!(decomp.location(), (2, 13));
        eprintln!("{:?}", decomp);
        assert_eq!(decomp.encoding().len(), 3);
    }
}
