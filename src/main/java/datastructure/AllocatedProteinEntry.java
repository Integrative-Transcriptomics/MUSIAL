package datastructure;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.stream.Collectors;

import components.Bio;

/**
 * Internal object representation of a protein that was allocated to a single {@link FeatureEntry}.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public class AllocatedProteinEntry {

    /**
     * {@link String}; The internal name of the protein.
     */
    public final String name;
    /**
     * {@link String} representation of a `pdb` file containing protein structure information of a coding gene feature.
     */
    public final String pdb;
    /**
     * Aligned reference amino acid sequences, one per protein chain, of the entry.
     */
    public final HashMap<String, String> chainSequences = new HashMap<>();
    /**
     * {@link String} identifier used to store the wild type proteoform.
     */
    public final static String WILD_TYPE_PROTEOFORM_ID = "WildType";
    /**
     * Map structure mapping {@link String}s (representing proteoform names) to {@link ProteoformEntry}s
     */
    public final ConcurrentSkipListMap<String, ProteoformEntry> proteoforms = new ConcurrentSkipListMap<>(String.CASE_INSENSITIVE_ORDER);
    /**
     * Hierarchical map structure that uses {@link String}s as keys and point to {@link AminoacidVariantEntry}s.
     * <p>
     * The first and second layer of the map structure use positions wrt. reference feature and the alternate amino-acid content as keys, respectively.
     */
    public final ConcurrentSkipListMap<String, ConcurrentSkipListMap<String, AminoacidVariantEntry>> variants =
            new ConcurrentSkipListMap<>((k1, k2) -> {
                int p1 = Integer.parseInt(k1.split("\\+")[0]);
                int p2 = Integer.parseInt(k2.split("\\+")[0]);
                if (p1 == p2) {
                    return Integer.compare(
                            k1.contains("+") ? k1.split("\\+")[1].hashCode() : 0,
                            k2.contains("+") ? k2.split("\\+")[1].hashCode() : 0
                    );
                }
                return Integer.compare(p1, p2);
            });

    /**
     * Constructor of {@link AllocatedProteinEntry}
     *
     * @param name           {@link String}; The internal name of the entry.
     * @param pdb            {@link String}; The file content of a .pdb format file.
     * @param chainSequences {@link HashMap}; Mapping chain ids (in the form of {@link String}s) to chain amino-acid sequences (in the form of {@link String}s).
     */
    public AllocatedProteinEntry(String name, String pdb, Map<String, String> chainSequences) {
        this.name = name;
        this.pdb = pdb;
        this.chainSequences.putAll(chainSequences);
    }

    /**
     * Returns a {@link String} that is used as internal proteoform name.
     * <p>
     * These names are generated by concatenating the prefix 'PF' with the HashCode of the {@link ProteoformEntry}
     * mandatory 'VSWAB' annotation value converted to a {@link String} for which a negative sign is converted to a one
     * and which in total was padded with prepending 0s to a length of 10.
     *
     * @param vSwab The {@link String} value of the 'VSWAB' annotation of a {@link ProteoformEntry}.
     * @return {@link String} intended to be used as internal proteoform name.
     */
    private String generateProteoformName(String vSwab) {
        if (vSwab.equals("")) {
            return AllocatedProteinEntry.WILD_TYPE_PROTEOFORM_ID;
        } else {
            StringBuilder proteoformNameBuilder = new StringBuilder();
            proteoformNameBuilder.append("PF");
            String hashCodeString = String.valueOf(vSwab.hashCode());
            if (hashCodeString.startsWith("-")) {
                proteoformNameBuilder.append("1");
                hashCodeString = hashCodeString.replace("-", "");
            } else {
                proteoformNameBuilder.append("0");
            }
            proteoformNameBuilder.append("0".repeat(10 - hashCodeString.length()));
            proteoformNameBuilder.append(hashCodeString);
            return proteoformNameBuilder.toString();
        }
    }

    /**
     * Adds a new {@link ProteoformEntry} to the {@link AllocatedProteinEntry#proteoforms} structure.
     *
     * @param fEntry   The {@link FeatureEntry} to whose {@link FeatureEntry#allocatedProtein} the entry shall be added.
     * @param sId      The {@link String} sample identifier of any {@link SampleEntry} of which the proteoform was derived.
     * @param variants {@link ConcurrentSkipListMap} of {@link String}s that maps positions wrt. the reference wild-type
     *                 proteoform to alternate amino-acid contents.
     */
    public void addProteoform(FeatureEntry fEntry, String sId, ConcurrentSkipListMap<String, String> variants) {
        String proteoformVSwab = variants.entrySet().stream().map(e -> e.getValue() + "@" + e.getKey()).collect(Collectors.joining("|"));
        float referenceProteinLength = (float) (fEntry.nucleotideSequence.length() / 3);
        DecimalFormat decimalFormat = new DecimalFormat("#.#");
        String proteoformName = generateProteoformName(proteoformVSwab);
        if (this.proteoforms.containsKey(proteoformName)) {
            this.proteoforms.get(proteoformName).samples.add(sId);
        } else {
            this.proteoforms.put(proteoformName, new ProteoformEntry(proteoformName, sId, proteoformVSwab));
            // FIXME: More efficient?
            // TODO: Extract as method.
            for (String variantPosition : variants.keySet()) {
                if (!this.variants.containsKey(variantPosition)) {
                    this.variants.put(variantPosition, new ConcurrentSkipListMap<>());
                }
                String variantContent = variants.get(variantPosition);
                // Check if variant induces an inner termination.
                if (variantContent.equals(String.valueOf(Bio.TERMINATION_AA1))) {
                    this.proteoforms.get(proteoformName).annotations.put("PT", "true");
                }
                if (!this.variants.get(variantPosition).containsKey(variantContent)) {
                    AminoacidVariantEntry aminoacidVariantAnnotationEntry = new AminoacidVariantEntry();
                    // TODO: Infer causative nucleotide variants.
                    this.variants.get(variantPosition).put(variantContent, aminoacidVariantAnnotationEntry);
                }
                this.variants.get(variantPosition).get(variantContent).occurrence.add(proteoformName);
            }
            if (!this.proteoforms.get(proteoformName).annotations.containsKey("PT")) {
                this.proteoforms.get(proteoformName).annotations.put("PT", "false");
            }
            // Compute percentage of variant positions; add as annotation VP.
            ArrayList<String> variantsContents = new ArrayList<>(variants.values());
            int variantPositions = variantsContents.subList(0, variantsContents.contains(String.valueOf(Bio.TERMINATION_AA1)) ? variantsContents.indexOf(String.valueOf(Bio.TERMINATION_AA1)) + 1 : variantsContents.size()).size();
            this.proteoforms.get(proteoformName).annotations.put("VP", decimalFormat.format(100 * (variantPositions / referenceProteinLength)).replace(",", "."));
        }

    }

}
