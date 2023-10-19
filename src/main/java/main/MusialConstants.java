package main;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Collection of common property keys used for annotations.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.2
 */
public class MusialConstants {
    public final static String REFERENCE_ID = "reference";
    public final static String VARIANTS = "variants";
    public final static String PROTEOFORM_DIFFERENTIAL_SEQUENCE = "differential_sequence";
    public final static String PROTEOFORM_NOVEL_STOPS = "novel_stops";
    public final static String NUMBER_SUBSTITUTIONS = "number_of_substitutions";
    public final static String NUMBER_INSERTIONS = "number_of_insertions";
    public final static String NUMBER_DELETIONS = "number_of_deletions";

    public final static String VARIABLE_POSITIONS = "variable_positions";

    public final static String FREQUENCY = "frequency";

    public final static String FREQUENCY_PASS = "frequency_pass";

    public final static String REFERENCE_CONTENT = "reference_content";

    public final static String SNP_EFF_PROPERTY_PREFIX = "snpeff_";
    public final static ArrayList<String> SNP_EFF_PROPERTIES = new ArrayList<>(Arrays.asList(
            "Effect",
            "Impact",
            "GeneName",
            "GeneID",
            "FeatureType",
            "FeatureID",
            "Biotype",
            "Rank/Total",
            "HGVS.c",
            "HGVS.p",
            "cDNA_position/cDNA_len",
            "CDS_position/CDS_len",
            "Protein_position/Protein_len"
    ));

    public final static String VARIANT_OCCURRENCE_SAMPLE_PREFIX = "of_sample_";

    public final static String SAMPLE_ANNOTATION_ALLELE_PREFIX = "allele_";

    public final static String SAMPLE_ANNOTATION_PROTEOFORM_PREFIX = "proteoform_";

    public final static String ALLELE_ANNOTATION_PROTEOFORM = "proteoform";

    public final static String FIELD_SEPARATOR_1 = ":";

    public final static String FIELD_SEPARATOR_2 = ";";
}
