package main;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Collection of common property keys used for annotations.
 *
 * @author Simon Hackl
 * @version 2.3
 * @since 2.2
 */
public class Constants {
    public final static String REFERENCE_FORM_NAME = "reference";
    public final static String NUMBER_OF_SUBSTITUTIONS = "number_of_substitutions";
    public final static String NUMBER_OF_INSERTIONS = "number_of_insertions";
    public final static String NUMBER_OF_DELETIONS = "number_of_deletions";
    public final static String NUMBER_OF_AMBIGUOUS = "number_of_ambiguous";
    public final static String VARIABLE_POSITIONS_PERCENTAGE = "variable_positions";
    public final static String FREQUENCY = "frequency";
    public final static String CONTENT_MODE_NUCLEOTIDE = "nucleotide";
    public final static String CONTENT_MODE_AMINOACID = "aminoacid";
    public final static String TYPE = "type";
    public final static String TYPE_AMBIGUOUS_PREFIX = "ambiguous_";
    public final static String TYPE_SUBSTITUTION = "substitution";
    public final static String TYPE_INSERTION = "insertion";
    public final static String TYPE_DELETION = "deletion";
    public final static String SNP_EFF_INFO_PREFIX = "snpeff_";
    public final static ArrayList<String> SNP_EFF_INFO_KEYS = new ArrayList<>(Arrays.asList(
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
    public final static String FIELD_SEPARATOR_1 = ":";
    public final static String FIELD_SEPARATOR_2 = ";";
    public final static String ANY_NUCLEOTIDE_STRING = "N";
    public final static char ANY_NUCLEOTiDE_CHAR = 'N';
    public final static String ANY_AMINOACID_STRING = "X";
    public final static String ANY_AMINOACID3_STRING = "ANY";
    public final static char ANY_AMINOACID_CHAR = 'X';
    public final static String TERMINATION_AMINOACID_STRING = "*";
    public final static String TERMINATION_AMINOACID3_STRING = "TER";
    public final static char TERMINATION_AMINOACID_CHAR = '*';
    public final static String DELETION_OR_GAP_STRING = "-";
    public final static String DELETION_OR_GAP_3_STRING = "GAP";
    public final static char DELETION_OR_GAP_CHAR = '-';
    public final static String CALL_INFO_NO_VARIANT = ".";
    public final static String CALL_INFO_NO_INFO = "?";
    public final static String CALL_INFO_REJECTED = "!";

}
