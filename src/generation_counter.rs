// Author: André Amado

//! Generation counter module

use crate::population::Population;

/// Struct that keeps track of generation and end conditions
pub struct GenerationCounter {
    t: usize,
    max_generation: Option<usize>
}

impl GenerationCounter {
    /// Creates a new generation counter
    ///
    /// # Arguments:
    /// * `max_time` 
    pub fn new(max_generation: Option<usize>) -> Self {
        Self { 
            t: 0,
            max_generation
        }
    }

    /// Advances the counter by one generation
    pub fn advance(&mut self) {
        self.t += 1;
    }

    /// Resets the counter to zero
    pub fn reset(&mut self) {
        self.t = 0;
    }

    /// Gets the current generation
    pub fn current(&self) -> usize {
        self.t
    }

    /// Checks if the run finished
    pub fn finished(&self, population: &Population, mutation: bool) -> bool {
        // If a maximum time is set check if it was reached
        if let Some(max_generation) = self.max_generation {
            if self.t >= max_generation {
                return true
            }
        }
        population.active_genotypes().len() == 1 && !mutation
    }

    /// Checks if maximum generation is set
    pub fn is_maximum_generation(&self) -> bool {
        self.max_generation.is_some()
    }
}