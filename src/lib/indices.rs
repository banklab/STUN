#![allow(dead_code)]

use serde::{Deserialize, Serialize};

use crate::lib::matrix::Matrix;

/// Struct that computes the convertions between indices and sequences of genotypes
#[derive(Serialize, Deserialize, Clone)]
pub struct Indices {
    l: usize,
    seqs: Matrix<bool>,
    powers: Vec<usize>,
}

impl Indices {
    /// Generates a new index structure for genotypes with `l` loci
    ///
    /// # Arguments:
    /// * `l` number of loci
    pub fn new(l: usize) -> Indices {
        let size = 2usize.pow(l as u32);
        let mut seqs = Matrix::<bool>::new(size, l);
        let mut powers = vec![0; l];

        for i in 0..l {
            powers[i] = 2usize.pow(i as u32);
        }

        for i in 0..size {
            let mut k = i;
            for j in 0..l {
                seqs.set(i, j, (k % 2) == 1);
                k = k / 2;
            }
        }

        let seqs = seqs;
        let powers = powers;

        Indices { l, seqs, powers }
    }

    /// Converts a sequence to an index
    ///
    /// # Arguments:
    /// * `seq` sequence for the genotype
    pub fn to_index(&self, seq: &[bool]) -> usize {
        let mut idx: usize = 0;
        for i in 0..self.l {
            idx += seq[i] as usize * self.powers[i];
        }
        idx
    }

    /// Converts an index to a sequence
    ///
    /// # Arguments:
    /// * `idx` index of the genotype
    pub fn to_sequence(&self, idx: usize) -> &[bool] {
        if idx > self.seqs.rows {
            panic!("index out of bounds {}. Maximum {}", idx, self.seqs.rows);
        }
        &self.seqs.get_row(idx)
    }

    /// Gets allele at locus `j` of genotype indexed `i`
    ///
    /// # Arguments:
    /// * `i` index of the genotype
    /// * `j` index of the locus
    pub fn get(&self, i: usize, j: usize) -> bool {
        self.seqs.get(i, j)
    }

    /// Gets an iterator that iterates over all genotype sequences
    pub fn seq_iter(&self) -> std::slice::Chunks<bool> {
        self.seqs.row_iter()
    }

    /// Gets the number of genotypes currently present in the population
    pub fn n_genotypes(&self) -> usize {
        self.seqs.rows
    }
}
