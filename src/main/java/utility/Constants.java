package utility;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * A utility class that holds constant values used throughout the project.
 * <p>
 * This class provides a centralized location for defining constants, ensuring
 * consistency and reducing the risk of hardcoding values in multiple places.
 * The constants include symbols, sequence types, and other project-specific
 * values. All fields are declared as {@code static final} to indicate that
 * they are constants and cannot be modified.
 */
public final class Constants {

    /**
     * Empty string constant.
     */
    public final static String EMPTY = "";

    /**
     * Colon string constant.
     */
    public final static String COLON = ":";

    /**
     * Semicolon string constant.
     */
    public final static String SEMICOLON = ";";

    /**
     * Comma string constant.
     */
    public final static String COMMA = ",";

    /**
     * Tabulator string constant.
     */
    public final static String TAB = "\t";

    /**
     * Dot string constant.
     */
    public final static String DOT = ".";

    /**
     * Equals sign string constant.
     */
    public final static String EQUAL = "=";

    /**
     * Number sign string constant.
     */
    public final static String SIGN = "#";

    /**
     * Pipe string constant.
     */
    public final static String PIPE = "|";

    /**
     * Represents the nucleotide used to denote any base in a sequence.
     * <p>
     * This constant is used in scenarios where a specific nucleotide cannot be determined
     * and is represented by the character 'N'.
     */
    public final static String anyNucleotide = "N";

    /**
     * Represents the stop codon in a sequence.
     * <p>
     * This constant is used to denote the stop codon, which is represented
     * by the character '*'. The stop codon signals the termination of
     * translation in a nucleotide sequence.
     */
    public final static String stopCodon = "*";

    /**
     * Prefix used to indicate low coverage variant calls.
     * <p>
     * This constant is used to mark variant calls with insufficient read depth,
     * represented by the prefix 'x'.
     */
    public final static String lowCoverageCallPrefix = "x";

    /**
     * Prefix used to indicate low frequency variant calls.
     * <p>
     * This constant is used to mark variant calls with low allele frequency,
     * represented by the prefix 'f'.
     */
    public final static String lowFrequencyCallPrefix = "f";

    /**
     * String representation of a gap in a sequence.
     */
    public final static String gapString = "-";

    /**
     * Character representation of a gap in a sequence.
     */
    public final static char gapChar = '-';

    /**
     * Represents the term used for synonymous {@link datastructure.Feature.Allele}s that are not actually stored.
     */
    public final static String synonymous = "synonymous";

    /**
     * Key used to represent the reference frequency of an attributable entity.
     * <p>
     * This constant is used as a key in data structures to store or retrieve
     * the reference frequency associated with a feature or sample.
     * <p>
     * For features, this corresponds to the proportion of samples that are
     * not associated with any variant on the feature locus. For samples,
     * this corresponds to the proportion of features for which the sample
     * is not associated with any variant.
     */
    public final static String $Attributable_frequencyReference = "frequency_reference";

    /**
     * Key used to represent the disrupted frequency of an attributable entity.
     * <p>
     * This constant is used as a key in data structures to store or retrieve
     * the reference frequency associated with a feature or sample and is
     * <b>only used if proteoform information is processed</b>.
     * <p>
     * For features, this corresponds to the proportion of proteoforms that are
     * associated with an {@code start_lost} or {@code stop_gained} SO effect.
     * For samples, this corresponds to the proportion of coding features for
     * which the sample is associated with a proteoform that is itself associated
     * with an {@code start_lost} or {@code stop_gained} SO effect.
     */
    public final static String $Attributable_frequencyDisrupted = "frequency_disrupted";

    /**
     * Key used to represent the length of a contig.
     */
    public final static String $Contig_length = "length";

    /**
     * Key used to represent the children of a feature.
     */
    public final static String $Feature_children = "children";

    /**
     * Key used to represent the number of alleles of a feature.
     * <p>
     * This constant is used as a key in data structures to store or retrieve
     * the number of alleles associated with a specific feature.
     */
    public final static String $Feature_numberOfAlleles = "no_allele";

    /**
     * Key used to represent the number of proteoforms of a feature.
     * <p>
     * This constant is used as a key in data structures to store or retrieve
     * the number of proteoforms associated with a specific feature.
     */
    public final static String $Feature_numberOfProteoforms = "no_proteoform";

    /**
     * Key used to represent the effects of a {@link datastructure.SequenceType}.
     * <p>
     * This constant is used as a key in
     * <ul>
     *     <li>{@link datastructure.Feature.Allele} data structures to store or retrieve
     *     the effects associated with a specific allele, i.e., all distinct {@code snpeff_effect} values of
     *     {@link datastructure.VariantInformation} instances associated with the object.</li>
     *     <li>{@link datastructure.Feature.Proteoform} data structures to store or retrieve
     *     the effects derived in the constructor of the class.</li>
     * </ul>
     */
    public final static String $SequenceType_effects = "so_effects";

    /**
     * Key used to represent the cluster label of a sequence type.
     */
    public final static String $SequenceType_clusterLabel = "cluster_label";

    /**
     * Key used to represent the cluster index of a sequence type.
     */
    public final static String $SequenceType_clusterIndex = "cluster_index";

    /**
     * Key used to represent the outlier score of a sequence type.
     */
    public final static String $SequenceType_clusterOutlierScore = "cluster_outlier_score";

    /**
     * Key used to represent the net shift in sequence length of a sequence type.
     * <p>
     * This constant is used as a key in data structures to store or retrieve
     * the net shift value associated with a specific sequence type. The net shift
     * typically represents the cumulative effect of insertions and deletions
     * on the sequence length.
     */
    public final static String $SequenceType_sequenceLengthVariation = "sequence_length_deviation";

    /**
     * Key used to represent the frequency of a sequence type.
     * <p>
     * This constant is used as a key in data structures to store or retrieve
     * the frequency value associated with a specific sequence type. The frequency
     * represents the proportion of samples that exhibit the specified sequence type.
     */
    public final static String $SequenceType_frequency = "allelic_frequency";

    /**
     * Key used to represent the proteoform of an allele.
     */
    public final static String $Allele_proteoform = "proteoform";

    /**
     * Key used to represent the number of calls in a sample.
     */
    public final static String $Sample_numberOfCalls = "no_call";

    /**
     * Key used to represent the number of filtered calls in a sample.
     */
    public final static String $Sample_numberOfFiltered = "no_call_filter";

    /**
     * Key used to represent the number of substitutions in a sample.
     */
    public final static String $Sample_numberOfSubstitutions = "no_variant_substitution";

    /**
     * Key used to represent the number of insertions in a sample.
     */
    public final static String $Sample_numberOfIndels = "no_variant_indel";

    /**
     * Key used to represent the mean coverage of a sample wrt. all calls.
     */
    public final static String $Sample_meanCoverage = "mean_coverage";

    /**
     * Key used to represent the mean coverage of a sample wrt. all calls.
     */
    public final static String $Sample_meanQuality = "mean_quality";

    /**
     * Key used to represent the frequency of a variant in the stored samples.
     */
    public final static String $VariantInformation_frequency = "variant_frequency";

    /**
     * The system-dependent line separator string.
     */
    public static final String lineSeparator = System.getProperty("line.separator");

    /**
     * The base symbols of the IUPAC nucleotide and amino acid code.
     */
    public static final String baseSymbols = "ARNDCQEGHILKMFPSTWYVBJZX*";

    /**
     * Prefix for SnpEff attribute keys.
     */
    public final static String snpEffAttributeKeyPrefix = "snpeff_";

    /**
     * SnpEff annotation field names.
     */
    public final static ArrayList<String> snpEffKeys = new ArrayList<>(Arrays.asList("allele", "effect", "impact", "gene_name", "gene_id", "feature_type", "feature_id", "biotype", "rank/total", "hgvs_c", "hgvs_p", "cDNA_position", "cds_position", "protein_position", "feature_distance", "note"));
}
