//! The `environment` module contains the setup of the evolutionary network.
pub use self::configuration::{Environment, EnvironmentBuilder};
pub use self::execution::GlobalEnvironment;

mod configuration;
mod execution;
