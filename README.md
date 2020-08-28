# oben (Organism Based Evolutionary Network)

A Rust framework for evolutionary approximation of arbitrary problems.

## Introduction

OBEN (german "oben": upside) stands for "organism based evolutionary network".
It simulates a population of organisms with varying genotypes and phenotypes.
Those organsims perform a user specified task, which is used to evalute their
fitness and thereby the amount of offspring they produce. Over time allele
combinations conveying a high fitness will increase in frequency, allowing
adaption to the specified task.

## How exactly the system works

A set of elementary information is defined, e.g. binary sequences.
This information represents the nodes of the network. Information can
flow through the network by elementary actions, termed reactions, for
example inversion of a binary sequence. Reactions are triggered by changing
information states, e.g. infomration A containing more set bits than B.
All aspects of this network can be subject to user defined mutation.

## Usage
