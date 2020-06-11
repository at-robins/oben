//! The `environment` module contains the setup of the evolutionary network.

use std::path::{Path, PathBuf};
use uuid::Uuid;

/// The sub-folder in which genome files are stored.
const SUBFOLDER_GENOME: &str = "genomes";
/// The sub-folder in which genome files of extinct populations are stored.
const SUBFOLDER_GENOME_EXTINCT: &str = "extinct";
/// The sub-folder in which population snapshot files are stored.
const SUBFOLDER_POPULATION: &str = "populations";
/// The file extension of genome files.
const FILE_EXTENSION_GENOME: &str = "genome";
/// The file extension of population files.
const FILE_EXTENSION_POPULATION: &str = "population";

/// An `Environment` specifing settings for an evolutionary network to develop in.
pub struct Environment {
    working_directory: PathBuf,
    /// The chance of a single offspring to carry a mutation.
    mutation_rate: f64,
    /// The maximum size a [`ClonalPopulation`] can grow to.
    ///
    /// [`ClonalPopulation`]: ./struct.ClonalPopulation.html
    max_clonal_population_size: u32,
    /// The rate in individualts / second at which individuals of a population die.
    death_rate: f64,
}

impl Environment {
    /// Returns the path to the working directory.
    pub fn working_directory(&self) -> &Path {
        Path::new(&self.working_directory)
    }

    /// Returns the chance of a single offspring carrying a mutation.
    pub fn mutation_rate(&self) -> f64 {
        self.mutation_rate
    }

    /// Returns the maximum size a [`ClonalPopulation`] can grow to.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    pub fn max_clonal_population_size(&self) -> u32 {
        self.max_clonal_population_size
    }

    /// Returns the rate in individualts / second at which individuals of a [`ClonalPopulation`] die.
    ///
    /// [`ClonalPopulation`]: ../population/struct.ClonalPopulation.html
    pub fn death_rate(&self) -> f64 {
        self.death_rate
    }

    /// Returns the file path to the [`Genome`] with the specified UUID.
    ///
    /// # Parameters
    ///
    /// * `genome_uuid` - the UUID of the [`Genome`]
    ///
    /// [`Genome`]: ../gene/struct.Genome.html
    pub fn genome_path(&self, genome_uuid: Uuid) -> PathBuf {
        let mut path_to_genome: PathBuf = self.working_directory().into();
        path_to_genome.push(SUBFOLDER_GENOME);
        path_to_genome.set_file_name(genome_uuid.to_string());
        path_to_genome.set_extension(FILE_EXTENSION_GENOME);
        path_to_genome
    }

}
