// author: André Amado

//! Matrix module

use serde::{Deserialize, Serialize};

/// Auxiliary type to store and manipulate simple matrices
#[derive(Serialize, Deserialize, Clone)]
pub struct Matrix<T> {
    data: Vec<T>,
    pub rows: usize,
    pub cols: usize,
}


impl<T: Default + Copy> Matrix<T> {
    /// Creates a new matrix
    ///
    /// # Arguments:
    /// * `rows` number of rows
    /// * `cols` number of columns
    pub fn new(rows: usize, cols: usize) -> Self {
        let mut data = Vec::<T>::new();
        data.resize_with(rows * cols, Default::default);
        Self { data, rows, cols }
    }

    /// Gets the ith row of the matrix
    ///
    /// # Arguments:
    /// * `i` index of the matrix row to retrieve
    #[inline]
    pub fn get_row(&self, i: usize) -> &[T] {
        &self.data[self.row(i)]
    }

    /// Gets the element (i, j) of the matrix
    ///
    /// # Arguments:
    /// * `i` row index
    /// * `j` column index
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> T {
        self.data[self.position(i, j)]
    }

    /// Sets the element (i, j) of the matrix to `element`
    ///
    /// # Arguments:
    /// * `i` row index
    /// * `j` column index
    /// * `element` value
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, element: T) {
        let position = self.position(i, j);
        self.data[position] = element;
    }

    /// Adds a row in the end of the matrix
    ///
    /// # Arguments:
    /// * `row` slice with the row contents to add
    pub fn add_row(&mut self, row: &[T]) {
        if row.len() != self.cols {
            panic!(
                "row length ({}) is not equal to the number of columns ({})!",
                row.len(),
                self.cols
            );
        }
        self.data.extend_from_slice(row);
        self.rows += 1;
    }

    /// Removes and returns the ith row of the matrix
    ///
    /// # Arguments:
    /// * `i` row index
    #[inline]
    pub fn extract_row(&mut self, i: usize) -> Vec<T> {
        self.data.drain(self.row(i)).collect()
    }

    /// Removes the ith row of the matrix replacing it with the last row of the matrix
    ///
    /// # Arguments:
    /// * `i` row index
    #[inline]
    pub fn swap_remove_row(&mut self, i: usize) {
        for j in self.row(i).rev() {
            self.data.swap_remove(j);
        }
    }

    /// Gets an iterator for the ith row of the matrix
    #[inline]
    pub fn row_iter(&self) -> std::slice::Chunks<T> {
        self.data.chunks(self.cols)
    }

    /// Gets the internal indices for the ith row of the matrix
    ///
    /// # Arguments:
    /// * `i` row index
    #[inline]
    fn row(&self, i: usize) -> std::ops::Range<usize> {
        (i * self.cols)..((i + 1) * self.cols)
    }

    /// Gets the internal index of element (i, j)
    ///
    /// # Arguments:
    /// * `i` row index
    /// * `j` column index
    #[inline]
    fn position(&self, i: usize, j: usize) -> usize {
        i * self.cols + j
    }
}
