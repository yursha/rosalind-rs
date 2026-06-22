/// Represents an unstructured pool of individual organisms categorized by their genotypes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndividualPopulation {
    pub homozygous_dominant: u64,
    pub heterozygous: u64,
    pub homozygous_recessive: u64,
}

impl IndividualPopulation {
    pub fn new(homozygous_dominant: u64, heterozygous: u64, homozygous_recessive: u64) -> Self {
        Self {
            homozygous_dominant,
            heterozygous,
            homozygous_recessive,
        }
    }

    /// Calculates the probability that two randomly selected organisms
    /// will produce offspring possessing at least one dominant allele.
    pub fn dominant_phenotype_probability(&self) -> f64 {
        let dom_homozygotes = self.homozygous_dominant as f64;
        let heterozygotes = self.heterozygous as f64;
        let rec_homozygotes = self.homozygous_recessive as f64;

        let total_population = dom_homozygotes + heterozygotes + rec_homozygotes;
        if total_population < 2.0 {
            return 0.0; // Avoid division by zero or mating a single organism
        }

        let num_ordered_pairs = total_population * (total_population - 1.0);

        // Calculate the probability of picking two parents that produce a homozygous recessive (aa) child:

        // 1. Both parents are homozygous recessive (aa + aa -> 100% aa offspring)
        let p_rec_from_rec_rec = (rec_homozygotes * (rec_homozygotes - 1.0)) / num_ordered_pairs;

        // 2. One parent is heterozygous, one is homozygous recessive (Aa + aa -> 50% aa offspring)
        // Note: this captures events "Aa + aa" and "aa + Aa" so we multiply it by 2
        let p_rec_from_het_rec = 2.0 * 0.5 * (heterozygotes * rec_homozygotes) / num_ordered_pairs;

        // 3. Both parents are heterozygous (Aa + Aa -> 25% aa offspring)
        let p_rec_from_het_het = 0.25 * (heterozygotes * (heterozygotes - 1.0)) / num_ordered_pairs;

        let p_rec = p_rec_from_rec_rec + p_rec_from_het_rec + p_rec_from_het_het;

        // Complement rule: P(at least one dominant) = 1 - P(both recessive)
        1.0 - p_rec
    }
}

/// Represents a population of explicit, structured breeding couples categorized by genotype pairings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CouplePopulation {
    pub dom_dom: u64, // AA-AA
    pub dom_het: u64, // AA-Aa
    pub dom_rec: u64, // AA-aa
    pub het_het: u64, // Aa-Aa
    pub het_rec: u64, // Aa-aa
    pub rec_rec: u64, // aa-aa
}

impl CouplePopulation {
    pub fn new(
        dom_dom: u64,
        dom_het: u64,
        dom_rec: u64,
        het_het: u64,
        het_rec: u64,
        rec_rec: u64,
    ) -> Self {
        Self {
            dom_dom,
            dom_het,
            dom_rec,
            het_het,
            het_rec,
            rec_rec,
        }
    }

    /// Calculates the expected number of offspring displaying the dominant phenotype
    /// under the assumption that every couple has exactly two offspring.
    pub fn expected_dominant_offspring(&self) -> f64 {
        (self.dom_dom as f64 * 2.0)
            + (self.dom_het as f64 * 2.0)
            + (self.dom_rec as f64 * 2.0)
            + (self.het_het as f64 * 1.5)
            + (self.het_rec as f64 * 1.0)
            + (self.rec_rec as f64 * 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper assertion function to handle floating-point precision tolerances
    fn assert_near(actual: f64, expected: f64) {
        let epsilon = 1e-5;
        assert!(
            (actual - expected).abs() < epsilon,
            "Expected {}, got {}",
            expected,
            actual
        );
    }

    #[test]
    fn test_individual_equal_distribution() {
        let population = IndividualPopulation::new(2, 2, 2);
        let probability = population.dominant_phenotype_probability();

        assert_near(probability, 0.78333);
    }

    #[test]
    fn test_individual_insufficient_population() {
        // Population of 1 cannot mate
        let population = IndividualPopulation::new(1, 0, 0);
        assert_near(population.dominant_phenotype_probability(), 0.0);

        // Empty population
        let empty_population = IndividualPopulation::new(0, 0, 0);
        assert_near(empty_population.dominant_phenotype_probability(), 0.0);
    }

    #[test]
    fn test_individual_pure_homozygous_dominant() {
        // If everyone is AA, probability of dominant offspring must be 1.0
        let population = IndividualPopulation::new(10, 0, 0);
        assert_near(population.dominant_phenotype_probability(), 1.0);
    }

    #[test]
    fn test_individual_pure_homozygous_recessive() {
        // If everyone is aa, probability of dominant offspring must be 0.0
        let population = IndividualPopulation::new(0, 0, 10);
        assert_near(population.dominant_phenotype_probability(), 0.0);
    }

    #[test]
    fn test_individual_dominant_and_recessive_only() {
        // AA (100) and aa (100).
        let population = IndividualPopulation::new(100, 0, 100);
        assert_near(population.dominant_phenotype_probability(), 0.75126);
    }

    #[test]
    fn test_couple_dominant_offspring_sample() {
        // Sample input: 1 0 0 1 0 1
        let couples = CouplePopulation::new(1, 0, 0, 1, 0, 1);
        assert_near(couples.expected_dominant_offspring(), 3.5);
    }

    #[test]
    fn test_couple_dominant_offspring_large() {
        let couples = CouplePopulation::new(18321, 19124, 16302, 15210, 17112, 19083);
        assert_near(couples.expected_dominant_offspring(), 147421.0);
    }
}
