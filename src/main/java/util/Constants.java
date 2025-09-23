package util;

import model.Allele;
import model.Proteoform;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Utility class that holds constant values used throughout the application.
 * <p>
 * This class provides a centralized location for defining constants, ensuring consistency and reducing the risk of hardcoding values in
 * multiple places. The constants include symbols, sequence types, and other project-specific values.
 */
public final class Constants {

    /**
     * Private constructor to prevent instantiation of this utility class.
     */
    private Constants() {
    }

    /**
     * Empty string constant.
     */
    public static final String EMPTY = "";

    /**
     * Colon string constant.
     */
    public static final String COLON = ":";

    /**
     * Semicolon string constant.
     */
    public static final String SEMICOLON = ";";

    /**
     * Comma string constant.
     */
    public static final String COMMA = ",";

    /**
     * Tabulator string constant.
     */
    public static final String TAB = "\t";

    /**
     * Dot string constant.
     */
    public static final String DOT = ".";

    /**
     * Equals sign string constant.
     */
    public static final String EQUAL = "=";

    /**
     * Number sign string constant.
     */
    public static final String SIGN = "#";

    /**
     * Underscore string constant.
     */
    public static final String UNDER_SCORE = "_";

    /**
     * Greater than sign string constant.
     */
    public static final String GREATER_THAN = ">";

    /**
     * Pipe symbol string constant.
     */
    public static final String PIPE = "\\|";

    /**
     * IUPAC symbol for any/unknown nucleotide in a nucleotide sequence.
     */
    public static final String ANY_NUCLEOTIDE = "N";

    /**
     * Translated stop codon in an amino acid sequence.
     */
    public static final String TERMINAL_AA = "*";

    /**
     * {@link String} representation of a gap in a sequence.
     */
    public static final String GAP = "-";

    /**
     * {@link Character} representation of a gap in a sequence.
     */
    public static final char GAP_CHAR = '-';

    /**
     * Prefix for genomic coordinate notation as per HGVS guidelines.
     */
    public static final String GENOMIC_COORDINATES_PREFIX = "g.";

    /**
     * Term used for reference alleles in the context of {@link Allele}s.
     */
    public static final String REFERENCE = "reference";

    /**
     * Term used for synonymous proteoforms in the context of {@link Proteoform}s.
     */
    public static final String SYNONYMOUS = "SYNONYMOUS";

    /**
     * Base symbols of the IUPAC nucleotide and amino acid code.
     */
    public static final String BASE_SYMBOLS = "ARNDCQEGHILKMFPSTWYVBJZX*";

    /**
     * Prefix for SnpEff attribute keys.
     */
    public static final String SNP_EFF_PREFIX = "snpeff_";

    /**
     * SnpEff annotation field names.
     * <p>
     * <ul>
     *     <li>0 = allele</li>
     *     <li>1 = effect</li>
     *     <li>2 = impact</li>
     *     <li>3 = gene_name</li>
     *     <li>4 = gene_id</li>
     *     <li>5 = feature_type</li>
     *     <li>6 = feature_id</li>
     *     <li>7 = biotype</li>
     *     <li>8 = rank/total</li>
     *     <li>9 = hgvs_c</li>
     *     <li>10 = hgvs_p</li>
     *     <li>11 = cDNA_position</li>
     *     <li>12 = cds_position</li>
     *     <li>13 = protein_position</li>
     *     <li>14 = feature_distance</li>
     *     <li>15 = note</li>
     * </ul>
     */
    public static final ArrayList<String> SNP_EFF_KEYS = new ArrayList<>(Arrays.asList("allele", "effect", "impact", "gene_name", "gene_id"
            , "feature_type", "feature_id", "biotype", "rank/total", "hgvs_c", "hgvs_p", "cDNA_position", "cds_position",
            "protein_position", "feature_distance", "note"));

    /**
     * System-dependent line separator string.
     */
    public static final String LINE_SEPARATOR = System.getProperty("line.separator");

    /**
     * A nested utility class that defines keys for various attributes used in the application.
     * <p>
     * This class provides a centralized location for defining attribute keys related to genomic data, ensuring consistency and reducing the
     * risk of hardcoding values in multiple places. These keys are used to store and retrieve metadata for samples, features, sequence
     * types, and variants.
     */
    public static final class AttributesKeys {

        /**
         * Key representing the fraction of reference alleles with respect to a {@link model.Sample} or {@link model.Feature}.
         */
        public static final String FREQUENCY_REFERENCE = "frequency_reference_allele";

        /**
         * Key representing the fraction of disrupted proteoforms with respect to a {@link model.Sample} or {@link model.Feature}.
         */
        public static final String FREQUENCY_DISRUPTED = "frequency_disrupted_proteoform";

        /**
         * Key representing the number of non-reference alleles associated with a feature.
         */
        public static final String NUMBER_OF_ALLELES = "count_allele";

        /**
         * Key representing the number of non-reference proteoforms associated with a feature.
         */
        public static final String NUMBER_OF_PROTEOFORMS = "count_proteoform";

        /**
         * Key representing the Sequence-Ontology effects associated with a {@link model.SequenceType}.
         */
        public static final String SO_EFFECTS = "so_effects";

        /**
         * Key representing the net shift in sequence length of a {@link model.SequenceType}.
         */
        public static final String SEQUENCE_LENGTH_DEVIATION = "length_delta";

        /**
         * Key representing the frequency of a {@link model.SequenceType} with respect to {@link model.Sample}s.
         */
        public static final String ALLELIC_FREQUENCY = "allelic_frequency";

        /**
         * Key representing the diversity of alleles associated with a {@link model.Feature}.
         */
        public static final String DIVERSITY_ALLELE = "diversity_allele";

        /**
         * Key representing the diversity of proteoforms associated with a {@link model.Feature}.
         */
        public static final String DIVERSITY_PROTEOFORM = "diversity_proteoform";

        /**
         * Key representing the fraction of filtered calls in a {@link model.Sample}.
         */
        public static final String FREQUENCY_FILTERED_CALLS = "frequency_calls_filtered";

        /**
         * Key representing the number of substitutions or single nucleotide variants (SNVs) in a {@link model.Sample}.
         */
        public static final String NUMBER_OF_SNVS = "no_variant_snv";

        /**
         * Key representing the number of insertions and deletions (InDels) in a {@link model.Sample}.
         */
        public static final String NUMBER_OF_INDELS = "no_variant_indel";

        /**
         * Key representing the mean coverage of a {@link model.Sample} over all of its calls.
         */
        public static final String MEAN_COVERAGE = "mean_coverage";

        /**
         * Key representing the mean entropy of a {@link model.Sample} over all of its calls.
         */
        public static final String MEAN_ENTROPY = "mean_entropy";

        /**
         * Key representing the frequency of a {@link model.Variant} across all stored samples.
         */
        public static final String VARIANT_FREQUENCY = "variant_frequency";
    }

}
