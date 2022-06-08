package datastructure;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.javatuples.Triplet;
import utility.Bio;

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
    public final static String WILD_TYPE_PROTEOFORM_ID = "WildType-0";
    public final ConcurrentSkipListMap<String, ProteoformEntry> proteoforms = new ConcurrentSkipListMap<>((k1, k2) -> {
        int i1 = Integer.parseInt(k1.split("-")[1]);
        int i2 = Integer.parseInt(k2.split("-")[1]);
        return Integer.compare(i1, i2);
    });
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
            return proteoformNameBuilder + "-" + vSwab.split("\\|").length;
        }
    }

    public void addProteoform(FeatureEntry fEntry, String sId,
                              ArrayList<Triplet<String, String, ArrayList<String>>> variableSegments) {
        String proteoformVSwab = variableSegments.stream().map(trplt -> trplt.getValue1() + "@" + trplt.getValue0()
        ).collect(Collectors.joining("|"));
        float referenceProteinLength = (float) (fEntry.nucleotideSequence.length() / 3);
        String proteoformName = generateProteoformName(proteoformVSwab);
        if (this.proteoforms.containsKey(proteoformName)) {
            this.proteoforms.get(proteoformName).samples.add(sId);
        } else {
            this.proteoforms.put(proteoformName, new ProteoformEntry(proteoformName, sId, proteoformVSwab));
            // FIXME: More efficient?
            for (Triplet<String, String, ArrayList<String>> variableSegment : variableSegments) {
                if (!this.variants.containsKey(variableSegment.getValue0())) {
                    this.variants.put(variableSegment.getValue0(), new ConcurrentSkipListMap<>());
                }
                int firstTerminationIndex = variableSegment.getValue1().indexOf(Bio.TERMINATION_AA1);
                if (firstTerminationIndex != -1) {
                    firstTerminationIndex += Integer.parseInt(variableSegment.getValue0().split("\\+")[0]);
                    if (firstTerminationIndex < referenceProteinLength) {
                        for (String chainSequence : fEntry.allocatedProtein.chainSequences.values()) {
                            if (chainSequence.substring(firstTerminationIndex - 1).chars()
                                    .anyMatch(c -> (char) c != Bio.TERMINATION_AA1 || (char) c != Bio.GAP ||
                                            (char) c != Bio.ANY_AA1)) {
                                this.proteoforms.get(proteoformName).annotations.put("PT", "true");
                                break;
                            }
                        }
                    }
                }
                if (!this.variants.get(variableSegment.getValue0()).containsKey(variableSegment.getValue1())) {
                    AminoacidVariantAnnotationEntry aminoacidVariantAnnotationEntry = new AminoacidVariantAnnotationEntry();
                    aminoacidVariantAnnotationEntry.annotations.put("cause",
                            String.join("|", variableSegment.getValue2()));
                    this.variants.get(variableSegment.getValue0()).put(variableSegment.getValue1(),
                            aminoacidVariantAnnotationEntry);
                }
                this.variants.get(variableSegment.getValue0()).get(variableSegment.getValue1()).occurrence.add(proteoformName);
            }
        }
    }

}
