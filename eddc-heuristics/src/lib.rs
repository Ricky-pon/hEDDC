use std::vec;

pub mod fasta;

pub struct NaiveEDDCParameters {
    gap_open_score: f64,
    gap_extention_score: f64,
    unit_duplication_score: Vec<f64>,
}

pub struct EDDCParameters {
    // Score of converting s[l, r) into u
    seg_conv_score: Vec<Vec<Vec<Vec<f64>>>>,
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

impl EDDCParameters {
    pub fn new(reads: &[Vec<usize>], units: &[&[u8]]) -> Self {
        use rayon::prelude::*;

        let unit_num = units.len();
        let unit_base_num = unit_num + 5; // units + ATGC + empty

        let mut unit_base = vec![];
        for &u in units {
            unit_base.push(u);
        }
        unit_base.push(b"A");
        unit_base.push(b"T");
        unit_base.push(b"G");
        unit_base.push(b"C");
        unit_base.push(b"");

        let naive_params = NaiveEDDCParameters::new(1.0, 1.0, units);

        let mut c1 = vec![vec![0.0; unit_base_num]; unit_base_num];
        for i in 0..unit_base_num {
            for j in i + 1..unit_base_num {
                let tmp = calc_c1(unit_base[i], unit_base[j], &unit_base, &naive_params);
                c1[i][j] = tmp;
                c1[j][i] = tmp;
            }
        }

        let mut c2 = vec![vec![vec![0.0; unit_base_num]; unit_base_num]; unit_base_num];
        for i in 0..unit_base_num {
            for j in 0..unit_base_num {
                for k in 0..unit_base_num {
                    let tmp = calc_c2(
                        unit_base[i],
                        unit_base[j],
                        unit_base[k],
                        &unit_base,
                        &naive_params,
                    );
                    c2[i][j][k] = tmp;
                }
            }
        }

        let mut valid_rule = vec![]; // (x, y, v): x y -> v
        let mut valid_rule_left = vec![]; // (x, v): x y -> v
        for x in 0..unit_base_num {
            for y in 0..unit_base_num {
                for v in 0..unit_base_num {
                    let mut flag = false;
                    for p in 0..unit_base_num {
                        for q in 0..unit_base_num {
                            if x == p && y == q {
                                continue;
                            }
                            if c2[x][y][v] >= c1[x][p] + c1[y][q] + c2[p][q][v] {
                                flag = true;
                                break;
                            }
                        }
                        if flag {
                            break;
                        }
                    }
                    if !flag {
                        valid_rule.push((x, y, v));
                        valid_rule_left.push((x, v));
                    }
                }
            }
        }
        valid_rule_left.sort();
        valid_rule_left.dedup();

        let seg_conv_score: Vec<Vec<Vec<Vec<f64>>>> = reads
            .par_iter()
            .map(|s| get_seg_conv_score(s, &c1, &c2, &valid_rule, &valid_rule_left))
            .collect();

        Self { seg_conv_score }
    }

    pub fn seg_conv_score(&self, seq_id: usize, l: usize, r: usize, unit_id: usize) -> f64 {
        self.seg_conv_score[seq_id][l][r][unit_id]
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

fn calc_c1(s: &[u8], t: &[u8], units: &[&[u8]], naive_params: &NaiveEDDCParameters) -> f64 {
    naive_eddc(s, t, units, naive_params)
}

fn calc_c2(
    s: &[u8],
    t: &[u8],
    u: &[u8],
    units: &[&[u8]],
    naive_params: &NaiveEDDCParameters,
) -> f64 {
    let mut st = vec![];
    for &c in s {
        st.push(c);
    }
    for &c in t {
        st.push(c);
    }

    naive_eddc(&st, u, units, naive_params)
}

fn get_seg_conv_score(
    s: &[usize],
    c1: &Vec<Vec<f64>>,
    c2: &[Vec<Vec<f64>>],
    valid_rule: &[(usize, usize, usize)],
    valid_rule_left: &[(usize, usize)],
) -> Vec<Vec<Vec<f64>>> {
    let n = s.len();
    let unit_num = c1.len();

    let mut dp = vec![vec![vec![f64::MAX; unit_num]; n + 1]; n + 1];
    let mut dp_sub = vec![vec![vec![vec![f64::MAX; unit_num]; unit_num]; n + 1]; n + 1];

    // initialization of dp
    for i in 0..n + 1 {
        for u in 0..unit_num {
            dp[i][i][u] = c1[unit_num - 1][u];

            if i < n {
                dp[i][i + 1][u] = c1[s[i]][u];
            }

            if i + 2 <= n {
                let mut val: Option<f64> = None;
                for a in 0..unit_num {
                    for b in 0..unit_num {
                        let tmp_val = c1[s[i]][a] + c1[s[i + 1]][b] + c2[a][b][u];
                        val = match val {
                            Some(v) => Some(v.min(tmp_val)),
                            None => Some(tmp_val),
                        };
                    }
                }
                dp[i][i + 2][u] = val.unwrap();
            }
        }
    }

    // initialization of dp_sub
    for i in 0..n + 1 {
        for j in i..(i + 3).min(n + 1) {
            for u in 0..unit_num {
                for a in 0..unit_num {
                    let mut val: Option<f64> = None;
                    for b in 0..unit_num {
                        let tmp_val = dp[i][j][b] + c2[a][b][u];
                        val = match val {
                            Some(v) => Some(v.min(tmp_val)),
                            None => Some(tmp_val),
                        };
                    }
                    dp_sub[i][j][u][a] = val.unwrap();
                }
            }
        }
    }

    // CYK-like algorithm
    for len in 3..n + 1 {
        for i in 0..n + 1 - len {
            let j = i + len;
            // dp
            for k in i + 1..j {
                for &(a, u) in valid_rule_left {
                    dp[i][j][u] = dp[i][j][u].min(dp[i][k][a] + dp_sub[k][j][u][a]);
                }
            }

            // dp_sub
            for &(a, b, u) in valid_rule {
                dp_sub[i][j][u][a] = dp_sub[i][j][u][a].min(dp[i][j][b] + c2[a][b][u]);
            }
        }
    }

    dp
}

fn alignment_dp_eddc(
    s_idx: usize,
    s: &[usize],
    t_idx: usize,
    t: &[usize],
    units: &[&[u8]],
    params: &EDDCParameters,
) -> f64 {
    let n = s.len();
    let m = t.len();
    let unit_num = units.len();
    let unit_base_num = unit_num + 5;
    let eps = unit_base_num - 1;

    let mut dp = vec![vec![0.0; m + 1]; n + 1];
    let mut dp_sub = vec![vec![vec![0.0; unit_base_num]; m + 1]; n + 1];

    // initialization of dp
    for j in 1..m + 1 {
        dp[0][j] = params.seg_conv_score(t_idx, 0, j, eps);
    }
    for i in 1..n + 1 {
        dp[i][0] = params.seg_conv_score(s_idx, 0, i, eps);
    }

    // initialization of dp_sub
    for j in 1..m + 1 {
        for u in 0..unit_base_num {
            let mut val: Option<f64> = None;
            for q in 0..j + 1 {
                let tmp_val = dp[0][q] + params.seg_conv_score(t_idx, q, j, u);
                val = match val {
                    Some(v) => Some(v.min(tmp_val)),
                    None => Some(tmp_val),
                };
            }
            dp_sub[0][j][u] = val.unwrap();
        }
    }

    // alignment
    for i in 1..n + 1 {
        for j in 1..m + 1 {
            // dp
            let mut val: Option<f64> = None;

            for u in 0..unit_base_num {
                // substitution & deletion
                for p in 0..i {
                    let tmp_val = params.seg_conv_score(s_idx, p, i, u) + dp_sub[p][j][u];
                    val = match val {
                        Some(v) => Some(v.min(tmp_val)),
                        None => Some(tmp_val),
                    };
                }
                // insertion
                for q in 0..j {
                    let tmp_val = dp[i][q]
                        + params.seg_conv_score(s_idx, i, i, u)
                        + params.seg_conv_score(t_idx, q, j, u);
                    val = match val {
                        Some(v) => Some(v.min(tmp_val)),
                        None => Some(tmp_val),
                    };
                }
            }

            dp[i][j] = val.unwrap();

            // dp_sub
            for u in 0..unit_base_num {
                let mut val: Option<f64> = None;

                for q in 0..j + 1 {
                    let tmp_val = dp[i][q] + params.seg_conv_score(t_idx, q, j, u);
                    val = match val {
                        Some(v) => Some(v.min(tmp_val)),
                        None => Some(tmp_val),
                    };
                }

                dp_sub[i][j][u] = val.unwrap();
            }
        }
    }

    dp[n][m]
}

pub fn eddc_heuristic(reads: &Vec<Vec<usize>>, units: &[&[u8]]) -> Vec<Vec<f64>> {
    let params = EDDCParameters::new(reads, units);

    let seq_num = reads.len();
    let mut scores = vec![vec![0.0; seq_num]; seq_num];

    for i in 0..seq_num {
        for j in i + 1..seq_num {
            let score = alignment_dp_eddc(i, &reads[i], j, &reads[j], units, &params);
            scores[i][j] = score;
            scores[j][i] = score;
        }
    }

    scores
}
/// Parallel version of the `eddc_heuristics`. The number of the thread is determined automatically.
/// To handle the # of the threads to be used, see [here](https://docs.rs/rayon/latest/rayon/struct.ThreadPoolBuilder.html#method.build_global).
pub fn eddc_heuristic_parallel(reads: &Vec<Vec<usize>>, units: &[&[u8]]) -> Vec<Vec<f64>> {
    let params = EDDCParameters::new(reads, units);

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
                    let score = alignment_dp_eddc(i, seq, j, seq_aligned, units, &params);
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn eddc_heuristic_test() {
        let sequences = vec![vec![0, 0, 0, 0, 0], vec![1, 1, 1, 1, 1]];
        let units = vec![b"AT".as_slice(), b"AAT".as_slice()];
        let correct_scores = vec![vec![0.0, 5.0], vec![5.0, 0.0]];
        let scores = eddc_heuristic(&sequences, &units);
        assert_eq!(correct_scores, scores);
    }
    #[test]
    fn eddc_heuristic_test2() {
        let sequences = vec![vec![0], vec![1, 1, 1]];
        let units = vec![b"ATTATTATC".as_slice(), b"ATT".as_slice()];
        let correct_scores = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let scores = eddc_heuristic(&sequences, &units);
        assert_eq!(correct_scores, scores);
    }
}
