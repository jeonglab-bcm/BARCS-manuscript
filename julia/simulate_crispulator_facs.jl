#!/usr/bin/env julia

# Simulate repeated CRISPR knockout FACS screens with Crispulator.jl 0.5.1.
#
# Each replicate starts from one shared guide library and produces:
#   low:    0--25% of the observed phenotype
#   high:   75--100%
#   bulk:   0--100%
#   input:  an independent sequencing draw from the bulk cell frequencies
#
# Low, bulk, and high are the requested three samples. Bulk overlaps both
# tails, so it is not an independent third FACS bin. Input and bulk deliberately
# carry no sorting direction and provide a negative-reference comparison.
#
# Optional positional arguments:
#   output_dir  replicates  seed  MOI  high_quality_guide_fraction  genes
# MOI, guide quality, and gene count can instead be supplied through
# CRISPULATOR_MOI, CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION, and
# CRISPULATOR_GENES.

using CSV
using Crispulator
using DataFrames
using DataStructures
using Distributions
using Random
using Statistics

const DEFAULT_OUTPUT = joinpath(
    dirname(@__DIR__), "data", "derived", "crispulator_facs"
)

function parse_integer(index::Int, default::Int)
    length(ARGS) < index && return default
    value = tryparse(Int, ARGS[index])
    isnothing(value) && error("Argument $index must be an integer")
    value
end

function parse_float(index::Int, default::Float64)
    length(ARGS) < index && return default
    value = tryparse(Float64, ARGS[index])
    isnothing(value) && error("Argument $index must be numeric")
    value
end

function environment_integer(name::String, default::Int)
    text = get(ENV, name, "")
    isempty(text) && return default
    value = tryparse(Int, text)
    isnothing(value) && error("Environment variable $name must be an integer")
    value
end

function environment_float(name::String, default::Float64)
    text = get(ENV, name, "")
    isempty(text) && return default
    value = tryparse(Float64, text)
    isnothing(value) && error("Environment variable $name must be numeric")
    value
end

output_dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : DEFAULT_OUTPUT
n_replicates = parse_integer(
    2, environment_integer("CRISPULATOR_REPLICATES", 4)
)
seed = parse_integer(3, 20250724)
moi = parse_float(4, environment_float("CRISPULATOR_MOI", 0.25))
high_quality_fraction = parse_float(
    5,
    environment_float(
        "CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION", 0.90
    ),
)
n_genes = parse_integer(
    6, environment_integer("CRISPULATOR_GENES", 400)
)

n_replicates >= 1 || error("At least one replicate is required")
0 < moi < 0.5 || error("MOI must be greater than zero and below 0.5")
0 <= high_quality_fraction <= 1 ||
    error("High-quality guide fraction must be between zero and one")
n_genes >= 20 || error("At least 20 genes are required")
mkpath(output_dir)

setup = FacsScreen()
setup.num_genes = n_genes
setup.coverage = environment_integer("CRISPULATOR_GUIDES_PER_GENE", 5)
setup.representation = environment_integer(
    "CRISPULATOR_TRANSFECTION_REPRESENTATION", 500
)
setup.moi = moi
setup.σ = environment_float("CRISPULATOR_PHENOTYPE_SD", 2.0)
setup.bottleneck_representation = environment_integer(
    "CRISPULATOR_SORT_REPRESENTATION", 50
)
setup.seq_depth = environment_integer("CRISPULATOR_READS_PER_GUIDE", 50)
setup.bin_info = OrderedDict(
    :low => (0.0, 0.25),
    :high => (0.75, 1.0),
    :bulk => (0.0, 1.0),
)

Random.seed!(seed)
phenotype_distributions =
    Dict{Symbol, Tuple{Float64, Distributions.Sampleable}}(
        :inactive => (0.75, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :increasing => (
            0.10, truncated(Normal(0.55, 0.20), 0.10, 1.0)
        ),
        :decreasing => (
            0.10, truncated(Normal(-0.55, 0.20), -1.0, -0.10)
        ),
    )
knockdown_distributions =
    Dict{Symbol, Tuple{Float64, Distributions.Sampleable}}(
        :high => (high_quality_fraction, Delta(1.0)),
        :low => (
            1 - high_quality_fraction,
            truncated(Normal(0.05, 0.07), 0.0, 1.0),
        ),
    )
library = if high_quality_fraction == 0.90
    # Preserve CRISPulator's exact default construction and random-number
    # sequence for the manuscript benchmark.
    Library(CRISPRn())
else
    Library(
        phenotype_distributions, knockdown_distributions, CRISPRn()
    )
end
guides, guide_frequency_distribution = construct_library(setup, library)
n_guides = length(guides)

guide_truth = DataFrame(
    guide = ["g$(lpad(index, 5, '0'))" for index in 1:n_guides],
    barcode_id = collect(1:n_guides),
    gene = ["GENE$(lpad(guide.gene, 4, '0'))" for guide in guides],
    gene_id = [guide.gene for guide in guides],
    class = string.([guide.class for guide in guides]),
    behavior = string.([guide.behavior for guide in guides]),
    knockdown = [guide.knockdown for guide in guides],
    theoretical_phenotype = [guide.theo_phenotype for guide in guides],
    library_frequency = guide_frequency_distribution.p,
)

count_table = DataFrames.select(guide_truth, [:guide, :gene])
sample_design = DataFrame(
    sample = String[],
    replicate = String[],
    sample_type = String[],
    lower_quantile = Union{Missing, Float64}[],
    upper_quantile = Union{Missing, Float64}[],
    directional_tail = Bool[],
)

sample_order = (:input, :low, :bulk, :high)
intervals = Dict(
    :input => (missing, missing),
    :low => (0.0, 0.25),
    :high => (0.75, 1.0),
    :bulk => (0.0, 1.0),
)

for replicate_index in 1:n_replicates
    # Separate seeds make every replicate independently reproducible while the
    # guide library and its ground-truth phenotypes remain shared.
    Random.seed!(seed + replicate_index)
    cells, cell_phenotypes = transfect(
        setup, library, guides, guide_frequency_distribution
    )
    sorted_cells = Crispulator.select(
        setup, cells, cell_phenotypes, guides
    )
    sample_frequencies = counts_to_freqs(sorted_cells, n_guides)

    # Input and bulk are independent sequencing draws from the same unsorted
    # transfected population. Their comparison must therefore be null apart
    # from sampling noise.
    sample_frequencies[:input] = copy(sample_frequencies[:bulk])
    depths = Dict(sample => setup.seq_depth for sample in sample_order)
    sequenced = sequencing(depths, guides, sample_frequencies)

    replicate_name = "R$(lpad(replicate_index, 2, '0'))"
    for sample_type in sample_order
        sample_name = "$(replicate_name)_$(sample_type)"
        count_table[!, sample_name] = sequenced[sample_type].counts
        interval = intervals[sample_type]
        push!(
            sample_design,
            (
                sample_name,
                replicate_name,
                string(sample_type),
                interval[1],
                interval[2],
                sample_type in (:low, :high),
            ),
        )
    end
end

gene_truth = combine(
    groupby(guide_truth, [:gene, :gene_id, :class, :behavior]),
    nrow => :n_guides,
    :theoretical_phenotype => median => :theoretical_phenotype,
)
gene_truth.active = in.(
    gene_truth.class, Ref(Set(["increasing", "decreasing"]))
)
gene_truth.expected_sign = ifelse.(
    gene_truth.class .== "increasing",
    1,
    ifelse.(gene_truth.class .== "decreasing", -1, 0),
)

parameters = DataFrame(
    parameter = [
        "crispulator_version",
        "seed",
        "replicates",
        "genes",
        "guides_per_gene",
        "cells_per_guide_at_transfection",
        "multiplicity_of_infection",
        "high_quality_guide_fraction",
        "phenotype_noise_sd",
        "sorted_cells_per_guide",
        "reads_per_guide_per_sample",
        "low_interval",
        "high_interval",
        "bulk_interval",
    ],
    value = string.([
        "0.5.1",
        seed,
        n_replicates,
        setup.num_genes,
        setup.coverage,
        setup.representation,
        setup.moi,
        high_quality_fraction,
        setup.σ,
        setup.bottleneck_representation,
        setup.seq_depth,
        "0-25%",
        "75-100%",
        "0-100%",
    ]),
)

CSV.write(joinpath(output_dir, "counts.tsv"), count_table; delim='\t')
CSV.write(joinpath(output_dir, "sample_design.tsv"), sample_design; delim='\t')
CSV.write(joinpath(output_dir, "guide_truth.tsv"), guide_truth; delim='\t')
CSV.write(joinpath(output_dir, "gene_truth.tsv"), gene_truth; delim='\t')
CSV.write(joinpath(output_dir, "parameters.tsv"), parameters; delim='\t')

println(
    "Wrote $(nrow(count_table)) guides, $(nrow(gene_truth)) genes, and ",
    "$(nrow(sample_design)) samples to $output_dir",
)
