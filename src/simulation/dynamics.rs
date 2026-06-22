/// Represents a generalized discrete-time population projection model.
/// Used to simulate population dynamics, epidemiological growth, and demographic shifts.
#[derive(Debug, Clone)]
pub struct AgeStructuredModel {
    /// The number of offspring produced per reproductive individual per time step.
    pub fecundity: u64,
    /// The probability of an individual surviving to the next time step (0.0 to 1.0).
    /// Set to 1.0 for immortal models; lower values simulate mortality.
    pub survival_rate: f64,
    /// If Some(m), cohorts at age m+1 are removed (senescence).
    pub max_lifespan: Option<usize>,
}

impl AgeStructuredModel {
    pub fn new(fecundity: u64, survival_rate: f64) -> Self {
        Self {
            fecundity,
            survival_rate: survival_rate.clamp(0.0, 1.0),
            max_lifespan: None,
        }
    }

    pub fn with_lifespan(mut self, lifespan: usize) -> Self {
        self.max_lifespan = Some(lifespan);
        self
    }

    /// Projects the population state forward by a given number of elapsed time intervals,
    /// starting from an initial population distribution.
    pub fn project(&self, starting_cohort_size: u128, elapsed_intervals: u32) -> u128 {
        // If mortal, we need `m` slots. If immortal, we only need 2 slots: [Newborns, Adults 2+]
        let cohort_count = self.max_lifespan.unwrap_or(2);

        let mut cohorts = vec![0u128; cohort_count];
        cohorts[0] = starting_cohort_size;

        // One interval is used up for the initial cohort to mature
        for _ in 1..elapsed_intervals {
            // 1. Offspring born from all mature cohorts (Index 1 to the end)
            let adults: u128 = cohorts[1..].iter().sum();
            let newborns = adults * (self.fecundity as u128);

            let mut next_cohorts = vec![0u128; cohort_count];

            // 2. Shift everyone forward one period
            for i in 0..(cohort_count - 1) {
                next_cohorts[i + 1] = self.survive(cohorts[i]);
            }

            // 3. THE MAGIC BIFURCATION: The Drain vs. The Cliff
            if self.max_lifespan.is_none() {
                // Immortal Mode: The elders in the last bucket survive AND stay in the last bucket
                next_cohorts[cohort_count - 1] += self.survive(cohorts[cohort_count - 1]);
            }
            // (If Some(m), we do nothing here. The oldest cohort simply gets overwritten/dropped!)

            next_cohorts[0] = newborns;
            cohorts = next_cohorts;
        }

        cohorts.iter().sum()
    }

    /// Helper to safely apply survival rates without ruining exact integer sequences
    fn survive(&self, count: u128) -> u128 {
        if (self.survival_rate - 1.0).abs() < f64::EPSILON {
            count // Bypass f64 entirely to protect exact Rosalind u128 precision
        } else {
            (count as f64 * self.survival_rate).round() as u128
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_fission_with_immortality() {
        // Standard binary fission / ideal growth (fecundity=1, survival=1.0)
        let model = AgeStructuredModel::new(1, 1.0);
        assert_eq!(model.project(10, 3), 20);
    }

    #[test]
    fn test_population_decay_with_mortality() {
        // No reproduction, 50% survival rate
        let model = AgeStructuredModel::new(0, 0.5);
        assert_eq!(model.project(1000, 1), 1000);
        assert_eq!(model.project(1000, 2), 500);
    }

    #[test]
    fn test_population_age_based_mortality() {
        let model = AgeStructuredModel::new(1, 1.0).with_lifespan(3);
        assert_eq!(model.project(1, 6), 4);
    }
}
