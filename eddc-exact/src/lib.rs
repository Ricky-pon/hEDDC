use std::vec;

pub mod fasta;

pub struct NaiveEDDCParameters {
    gap_open_score: f64,
    gap_extention_score: f64,
    unit_duplication_score: Vec<f64>,
}

impl NaiveEDDCParameters {
    pub fn new(gap_open_score: f64, gap_extention_score: f64, units: &[&[u8]]) -> Self {
        let unit_duplication_score: Vec<f64> =
            units.iter().map(|u| unit_duplication_score(u)).collect();

        Self {
            gap_open_score,
            gap_extention_score,
            unit_duplication_score,
        }
    }

    pub fn gap_open_score(&self) -> f64 {
        self.gap_open_score
    }

    pub fn gap_extention_score(&self) -> f64 {
        self.gap_extention_score
    }

    pub fn unit_duplication_score(&self, a: usize) -> f64 {
        self.unit_duplication_score[a]
    }
}

fn edit_distance(u: &[u8], v: &[u8]) -> f64 {
    let n = u.len();
    let m = v.len();
    let mut dp: Vec<Vec<Option<f64>>> = vec![vec![None; m + 1]; n + 1];

    dp[0][0] = Some(0.0);
    for i in 0..n + 1 {
        for j in 0..m + 1 {
            if i < n {
                let tmp_val = dp[i][j].unwrap() + 1.0;
                dp[i + 1][j] = match dp[i + 1][j] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if j < m {
                let tmp_val = dp[i][j].unwrap() + 1.0;
                dp[i][j + 1] = match dp[i][j + 1] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            if i < n && j < m {
                let tmp_val = dp[i][j].unwrap() + (if u[i] != v[j] { 1.0 } else { 0.0 });
                dp[i + 1][j + 1] = match dp[i + 1][j + 1] {
                    Some(val) => Some(val.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
        }
    }

    dp[n][m].unwrap()
}

fn unit_duplication_score(unit: &[u8]) -> f64 {
    unit.len() as f64 * 0.5
}

fn calc_s_to_a(
    s: &[u8],
    units: &[&[u8]],
    params: &NaiveEDDCParameters,
) -> Vec<Vec<Vec<Option<f64>>>> {
    let n = s.len();
    let u_num = units.len() - 5;
    let mut s_to_a: Vec<Vec<Vec<Option<f64>>>> = vec![vec![vec![None; n + 1]; n + 1]; u_num + 5];

    // [i, j) => unit
    for a in 0..u_num {
        for i in 0..n + 1 {
            s_to_a[a][i][i] = Some(params.unit_duplication_score(a));
        }
    }
    for a in 0..u_num {
        for len in 1..n + 1 {
            for i in 0..(n + 1 - len) {
                let j = i + len;
                let mut val = edit_distance(&s[i..j], units[a]);
                for k in i + 1..j {
                    val = val.min(
                        params.unit_duplication_score(a)
                            + s_to_a[a][i][k].unwrap()
                            + s_to_a[a][k][j].unwrap(),
                    );
                }
                s_to_a[a][i][j] = Some(val);
            }
        }
    }
    // [i, j) -> eps
    for i in 0..n + 1 {
        s_to_a[u_num + 4][i][i] = Some(0.0);
    }
    for len in 1..n + 1 {
        for i in 0..(n + 1 - len) {
            let j = i + len;
            let mut val = params.gap_open_score + params.gap_extention_score * (len - 1) as f64;
            for a in 0..u_num {
                val = val.min(s_to_a[a][i][j].unwrap() + params.unit_duplication_score(a));
            }
            for k in i + 1..j {
                val = val.min(s_to_a[u_num + 4][i][k].unwrap() + s_to_a[u_num + 4][k][j].unwrap());
            }
            s_to_a[u_num + 4][i][j] = Some(val);
        }
    }
    // [i, j) -> base
    for a in u_num..u_num + 4 {
        for i in 0..n + 1 {
            s_to_a[a][i][i] = Some(params.gap_open_score());
        }
        for i in 0..n {
            s_to_a[a][i][i + 1] = Some(if s[i] == units[a][0] { 0.0 } else { 1.0 });
        }
    }
    for a in u_num..u_num + 4 {
        for len in 2..n + 1 {
            for i in 0..(n + 1 - len) {
                let j = i + len;
                let mut val: Option<f64> = None;
                for k in i + 1..j {
                    let mut tmp_val = s_to_a[a][i][k].unwrap() + s_to_a[u_num + 4][k][j].unwrap();
                    tmp_val =
                        tmp_val.min(s_to_a[u_num + 4][i][k].unwrap() + s_to_a[a][k][j].unwrap());
                    val = match val {
                        Some(v) => Some(v.min(tmp_val)),
                        None => Some(tmp_val),
                    }
                }
                s_to_a[a][i][j] = val;
            }
        }
    }

    s_to_a
}

fn naive_eddc(s: &[u8], t: &[u8], units: &[&[u8]], params: &NaiveEDDCParameters) -> f64 {
    // units: U + Sigma + epsilon
    let n = s.len();
    let m = t.len();
    let u_num = units.len() - 5;

    let s_to_a = calc_s_to_a(s, units, params);
    let t_to_a = calc_s_to_a(t, units, params);

    let mut ed: Vec<Vec<Option<f64>>> = vec![vec![None; m + 1]; n + 1];
    let mut ed2: Vec<Vec<Vec<Option<f64>>>> = vec![vec![vec![None; m + 1]; n + 1]; u_num + 5];

    // initialization
    for j in 0..m + 1 {
        ed[0][j] = t_to_a[u_num + 4][0][j];
        for a in 0..u_num + 5 {
            let mut val: Option<f64> = None;
            for q in 0..j {
                let tmp_val = ed[0][q].unwrap() + t_to_a[a][q][j].unwrap();
                val = match val {
                    Some(v) => Some(v.min(tmp_val)),
                    None => Some(tmp_val),
                }
            }
            ed2[a][0][j] = val;
        }
    }
    for i in 0..n + 1 {
        ed[i][0] = s_to_a[u_num + 4][0][i];
    }
    // dp
    for i in 1..n + 1 {
        // ed
        for j in 1..m + 1 {
            let mut val: Option<f64> = None;
            for a in 0..u_num + 5 {
                for p in 0..i {
                    let tmp_val = s_to_a[a][p][i].unwrap() + ed2[a][p][j].unwrap();
                    val = match val {
                        Some(v) => Some(v.min(tmp_val)),
                        None => Some(tmp_val),
                    }
                }
                for q in 0..j {
                    let tmp_val =
                        ed[i][q].unwrap() + s_to_a[a][i][i].unwrap() + t_to_a[a][q][j].unwrap();
                    val = match val {
                        Some(v) => Some(v.min(tmp_val)),
                        None => Some(tmp_val),
                    }
                }
            }
            ed[i][j] = val;
        }

        // ed2
        for j in 1..m + 1 {
            for a in 0..u_num + 5 {
                let mut val: Option<f64> = None;
                for q in 0..j + 1 {
                    let tmp_val = ed[i][q].unwrap() + t_to_a[a][q][j].unwrap();
                    val = match val {
                        Some(v) => Some(v.min(tmp_val)),
                        None => Some(tmp_val),
                    }
                }
                ed2[a][i][j] = val;
            }
        }
    }

    ed[n][m].unwrap()
}

pub fn eddc_exact_parallel(reads: &[&[u8]], units: &[&[u8]]) -> Vec<Vec<f64>> {
    let mut unit_base = vec![];
    for &u in units {
        unit_base.push(u);
    }
    unit_base.push(b"A");
    unit_base.push(b"T");
    unit_base.push(b"G");
    unit_base.push(b"C");
    unit_base.push(b"");

    let params = NaiveEDDCParameters::new(1.0, 1.0, units);

    use rayon::prelude::*;
    let seq_num = reads.len();
    let scores: Vec<_> = reads
        .par_iter()
        .enumerate()
        .flat_map(|(i, seq)| {
            reads
                .iter()
                .enumerate()
                .skip(i + 1)
                .map(|(j, seq_aligned)| {
                    let score = naive_eddc(seq, seq_aligned, &unit_base, &params);
                    (i, j, score)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let mut score_matrix = vec![vec![0.0; seq_num]; seq_num];
    for (i, j, score) in scores {
        score_matrix[i][j] = score;
        score_matrix[j][i] = score;
    }
    score_matrix
}
