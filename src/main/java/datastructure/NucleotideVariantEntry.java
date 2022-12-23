package datastructure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ConcurrentSkipListSet;

/**
 * Object used to store annotations for a single nucleotide variant.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class NucleotideVariantEntry {

    /**
     * {@link HashMap} of {@link String} key/value pairs specifying any information related to the variant.
     */
    public final HashMap<String, String> annotations = new HashMap<>();

    /**
     * {@link ConcurrentSkipListSet} of {@link String}s indicating the occurrence of this variant.
     */
    public final ConcurrentSkipListSet<String> occurrence = new ConcurrentSkipListSet<>();
    public final static String PROPERTY_NAME_PRIMARY = "PRIMARY";
    public final static String PROPERTY_NAME_REFERENCE_CONTENT = "REF_CONT";

    public final static String PROPERTY_PREFIX_SNP_EF = "SNPEFF_";
    public final static ArrayList<String> PROPERTY_LIST_SNP_EFF = new ArrayList<>(Arrays.asList(
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
    public final static String PROPERTY_NAME_FREQUENCY = "FRQ";
    public final static String PROPERTY_NAME_LOCATION = "LOC";

    /**
     * Constructor of {@link NucleotideVariantEntry}.
     */
    public NucleotideVariantEntry() {

    }

}