// fn main() {
//     println!(
//         "{}",
//         compute_mcclusky_bitopt(
//             vec![
//                 "1".into(),
//                 "2".into(),
//                 "3".into(),
//                 "4".into(),
//                 "5".into(),
//                 "6".into(),
//                 "7".into(),
//                 "8".into(),
//                 "9".into(),
//                 "10".into(),
//                 "11".into(),
//                 "12".into(),
//                 "13".into(),
//                 "14".into(),
//             ],
//             vec![0, 1],
//             (4..2u32.pow(14)).collect(),
//         )
//     );
// }

use rayon::prelude::*;
use std::{
    collections::HashSet,
    time::{Duration, Instant},
};

type Bits = (u32, u32);

#[derive(Debug, Clone)]
struct MintermTableRow {
    minterms: Vec<u32>,
    bits: Bits,
    combined: bool,
}

impl MintermTableRow {
    fn high_bits_count(&self) -> usize {
        self.bits.0.count_ones() as usize
    }
}

type MintermTable = Vec<MintermTableRow>;

pub fn compute_mcclusky(
    input_names: Vec<String>,
    minterms: Vec<u32>,
    dont_cares: Vec<u32>,
) -> String {
    if minterms.is_empty() {
        return "0".to_string();
    }

    if minterms.len() + dont_cares.len() == 2usize.pow(input_names.len() as u32) {
        return "1".to_string();
    }

    let mut get_pairs_time: Duration = Duration::ZERO;
    let mut apply_pairs_time: Duration = Duration::ZERO;

    // Step 1: finding prime implicants

    // create first minterm table
    let mut minterm_tables: Vec<MintermTable> = vec![minterms
        .iter()
        .chain(dont_cares.iter())
        .map(|&mt| MintermTableRow {
            minterms: vec![mt],
            combined: false,
            bits: (mt, 0),
        })
        .collect()];
    //TODO: impl mem swap
    // let mut mt1: Vec<MintermTableRow> = vec![];
    // let mut mt2: Vec<MintermTableRow> = vec![];
    // std::mem::swap(&mut mt1, &mut mt2);

    let mut iteration: usize = 1;
    loop {
        println!(
            "Iteration: {}. Previous iteration produced {} rows.",
            iteration,
            minterm_tables[iteration - 1].len()
        );

        // create new table
        minterm_tables.push(vec![]);

        // sort table based on number of ones
        minterm_tables[iteration - 1].sort_by_key(|a| a.high_bits_count());

        // create group start indexes
        let group_start_indexes: Vec<Option<usize>> = (0..input_names.len() + 1)
            .map(|g| {
                minterm_tables[iteration - 1]
                    .iter()
                    .position(|row| row.high_bits_count() == g)
            })
            .collect();

        // iterate over groups
        for group in 0..group_start_indexes.len() - 1 {
            let Some(first_group_start_index) = group_start_indexes[group] else {
                continue;
            };

            let Some(second_group_start_index) = group_start_indexes[group + 1] else {
                continue;
            };

            let second_group_end_index = group_start_indexes
                .iter()
                .skip(group + 2)
                .find_map(|&gi| gi)
                .unwrap_or(minterm_tables[iteration - 1].len());

            let first_group_slice = first_group_start_index..second_group_start_index;
            let second_group_slice = second_group_start_index..second_group_end_index;

            // let mut pairs_to_combine: Vec<(usize, usize, usize)> = Vec::new();
            // for (i, row) in minterm_tables[iteration - 1][first_group_slice]
            //     .iter()
            //     .enumerate()
            // {
            //     for (j, other_row) in minterm_tables[iteration - 1][second_group_slice.clone()]
            //         .iter()
            //         .enumerate()
            //     {
            //         // check if bit difference is one
            //         let bit_differences = bit_differences(&row.bits, &other_row.bits);
            //         if bit_differences.len() == 1 {
            //             pairs_to_combine.push((i, j, bit_differences[0] as usize));
            //         }
            //     }
            // }
            let mut before = Instant::now();
            let pairs_to_combine: Vec<(usize, usize, u32)> = minterm_tables[iteration - 1]
                [first_group_slice]
                .par_iter()
                .enumerate()
                .flat_map(|(i, row)| {
                    minterm_tables[iteration - 1][second_group_slice.clone()]
                        .iter()
                        .enumerate()
                        .filter_map(|(j, other_row)| {
                            // bit_diff(&row.bits, &other_row.bits).map(|diff| (i, j, diff as usize))
                            if row.bits.1 & other_row.bits.1 != row.bits.1 {
                                return None;
                            }
                            let diff = row.bits.0 ^ other_row.bits.0;
                            if diff.count_ones() == 1 {
                                Some((i, j, diff.ilog2()))
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<(usize, usize, u32)>>()
                })
                .collect();

            get_pairs_time += before.elapsed();
            before = Instant::now();

            for (i, j, diff_index) in pairs_to_combine {
                let first_row = &minterm_tables[iteration - 1][first_group_start_index + i];
                let second_row = &minterm_tables[iteration - 1][second_group_start_index + j];

                let combined_minterms: Vec<u32> = first_row
                    .minterms
                    .iter()
                    .chain(
                        second_row
                            .minterms
                            .iter()
                            .filter(|mt| !first_row.minterms.contains(mt)),
                    )
                    .cloned()
                    .collect();

                let mut combined_bits: Bits = first_row.bits;
                combined_bits.1 ^= 2u32.pow(diff_index);

                // create combined row
                minterm_tables[iteration].push(MintermTableRow {
                    minterms: combined_minterms,
                    bits: combined_bits,
                    combined: false,
                });

                // mark both rows as combined
                minterm_tables[iteration - 1][first_group_start_index + i].combined = true;
                minterm_tables[iteration - 1][second_group_start_index + j].combined = true;
            }

            apply_pairs_time += before.elapsed();
        }

        // remove duplicate prime implicants
        let mut seen: HashSet<Bits> = HashSet::new();
        minterm_tables[iteration].retain(|row| {
            if seen.contains(&row.bits) {
                false
            } else {
                seen.insert(row.bits);
                true
            }
        });

        // check if finished combining
        if minterm_tables[iteration].is_empty() {
            break;
        }

        iteration += 1;
    }

    // Step 2: prime implicant chart
    let prime_implicants: Vec<MintermTableRow> = minterm_tables
        .iter()
        .flatten()
        .filter(|row| !row.combined)
        .cloned()
        .collect();

    let mut essential_prime_implicants: Vec<MintermTableRow> = prime_implicants
        .iter()
        .filter(|pi| {
            pi.minterms.iter().any(|mt| {
                !dont_cares.contains(mt)
                    && prime_implicants
                        .iter()
                        .flat_map(|pi2| &pi2.minterms)
                        .filter(|&mt2| mt2 == mt)
                        .count()
                        == 1
            })
        })
        .cloned()
        .collect();

    let not_covered_minterms: Vec<u32> = minterms
        .iter()
        .filter(|mt| {
            !essential_prime_implicants
                .iter()
                .flat_map(|pi| &pi.minterms)
                .any(|mt2| &mt2 == mt)
        })
        .cloned()
        .collect();

    if !not_covered_minterms.is_empty() {
        let prime_implicants_power_set = powerset(&prime_implicants);

        let mut covering_prime_implicant_sets: Vec<Vec<MintermTableRow>> =
            prime_implicants_power_set
                .iter()
                .filter(|set| {
                    let covered_minterms: Vec<u32> =
                        set.iter().flat_map(|row| &row.minterms).cloned().collect();

                    not_covered_minterms
                        .iter()
                        .all(|mt| covered_minterms.contains(mt))
                })
                .cloned()
                .collect();

        // sort by length of set first and then by most amount of minterms
        covering_prime_implicant_sets.sort_by_key(|s| s.len());
        let shortest_set_length = covering_prime_implicant_sets[0].len();
        covering_prime_implicant_sets.retain(|set| set.len() == shortest_set_length);
        covering_prime_implicant_sets
            .sort_by_key(|s| s.iter().flat_map(|row| &row.minterms).count());
        covering_prime_implicant_sets.reverse();

        essential_prime_implicants.append(&mut covering_prime_implicant_sets[0]);
    }

    println!("Get Pairs took: {:.2?}", get_pairs_time);
    println!("Apply Pairs took: {:.2?}", apply_pairs_time);

    to_expression_string(&essential_prime_implicants, input_names)
}

fn to_expression_string(rows: &[MintermTableRow], input_names: Vec<String>) -> String {
    rows.iter()
        .map(|row| {
            let conjunction: String = (0u32..input_names.len() as u32)
                .rev()
                .filter_map(|bit_index| {
                    if (row.bits.1 & 2u32.pow(bit_index)) > 0 {
                        None
                    } else if row.bits.0 & 2u32.pow(bit_index) > 0 {
                        Some(input_names[(input_names.len() - 1) - bit_index as usize].clone())
                    } else {
                        Some(format!(
                            "¬{}",
                            input_names[(input_names.len() - 1) - bit_index as usize]
                        ))
                    }
                })
                .collect::<Vec<String>>()
                .join(" ∧ ");

            format!("({})", conjunction)
        })
        .collect::<Vec<String>>()
        .join(" ∨ ")
}

fn powerset<T>(s: &[T]) -> Vec<Vec<T>>
where
    T: Clone,
{
    (0..2usize.pow(s.len() as u32))
        .map(|i| {
            s.iter()
                .enumerate()
                .filter(|&(t, _)| (i >> t) % 2 == 1)
                .map(|(_, element)| element.clone())
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::compute_mcclusky;

    #[test]
    fn trivial_zero() {
        assert_eq!(
            compute_mcclusky(vec!["S0".into(), "S1".into()], vec![], vec![0, 1]),
            "0"
        );
    }

    #[test]
    fn trivial_one() {
        assert_eq!(
            compute_mcclusky(vec!["S0".into(), "S1".into()], vec![0, 1, 2], vec![3]),
            "1"
        );
    }

    #[test]
    fn three_vars() {
        assert_eq!(
            compute_mcclusky(
                vec!["S0".into(), "S1".into(), "S2".into()],
                vec![0, 3, 4],
                vec![2, 7]
            ),
            "(¬S1 ∧ ¬S2) ∨ (S1 ∧ S2)"
        );
    }

    #[test]
    fn four_vars() {
        assert_eq!(
            compute_mcclusky(
                vec!["S0".into(), "S1".into(), "S2".into(), "S3".into()],
                vec![0, 3, 4, 8, 9, 10, 12, 15],
                vec![2, 7, 11]
            ),
            "(¬S2 ∧ ¬S3) ∨ (S0 ∧ ¬S1) ∨ (S2 ∧ S3)"
        );
    }
}
