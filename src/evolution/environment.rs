//! The `environment` module contains the setup of the evolutionary network.
pub use self::configuration::{Environment, EnvironmentBuilder};
pub use self::execution::EcologicalNiche;
pub use self::mutation::{Mutation, MutationCompendium};

mod configuration;
mod execution;
mod mutation;
