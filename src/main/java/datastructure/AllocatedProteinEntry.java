package datastructure;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.stream.Collectors;

import components.Bio;

/**
 * TODO
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public class AllocatedProteinEntry {

    /**
     *
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
    public final ConcurrentSkipListMap<String, ProteoformEntry> proteoforms = new ConcurrentSkipListMap<>(String.CASE_INSENSITIVE_ORDER);
    /**
     *
     */
    public final ConcurrentSkipListMap<String, ConcurrentSkipListMap<String, AminoacidVariantAnnotationEntry>> variants =
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

    public AllocatedProteinEntry(String name, String pdb, Map<String, String> chainSequences) {
        this.name = name;
        this.pdb = pdb;
        this.chainSequences.putAll(chainSequences);
    }

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
                    AminoacidVariantAnnotationEntry aminoacidVariantAnnotationEntry = new AminoacidVariantAnnotationEntry();
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
