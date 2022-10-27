package datastructure;

import components.Bio;
import components.Logging;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;

import java.util.*;
import java.util.stream.Collectors;

import static datastructure.VariantsDictionary.FIELD_SEPARATOR_1;

public class VariantsTable {

    private TreeMap<String, HashMap<String, String>> variantsTable = new TreeMap<>((s1, s2) -> {
        int p1 = Integer.parseInt(s1.split("\\+")[0]);
        int p2 = Integer.parseInt(s2.split("\\+")[0]);
        if (p1 != p2) {
            return Integer.compare(p1, p2);
        } else {
            int i1 = s1.contains("\\+") ? Integer.parseInt(s1.split("\\+")[1]) : 0;
            int i2 = s2.contains("\\+") ? Integer.parseInt(s2.split("\\+")[1]) : 0;
            return Integer.compare(i1, i2);
        }
    });

    private VariantsDictionary parentDictionary;

    private ArrayList<String> sampleIdentifiers;

    private String featureIdentfier;

    public final String mode;

    public final boolean excludeIndels;

    public final boolean includeNonVariantPositions;

    public final boolean redundantReferenceContent;

    public static final String MODE_NUCLEOTIDE = "NUCLEOTIDE";

    public static final String MODE_AMINOACID = "AMINOACID";

    public static final String LAYER_SAMPLE = "SAMPLE";

    public static final String LAYER_GROUPED = "GROUPED";

    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(Variants Table Construction Failed)";

    public VariantsTable(VariantsDictionary parentDictionary, ArrayList<String> sampleIdentifiers,
                         String featureIdentfier, String mode, boolean excludeIndels,
                         boolean includeNonVariantPositions, boolean redundantReferenceContent) throws MusialException {
        this.parentDictionary = parentDictionary;
        this.sampleIdentifiers = this.parentDictionary.removeInvalidSampleIds(sampleIdentifiers);
        this.featureIdentfier = featureIdentfier;
        this.mode = mode;
        this.excludeIndels = excludeIndels;
        this.includeNonVariantPositions = includeNonVariantPositions;
        this.redundantReferenceContent = redundantReferenceContent;
        // Check validity of specified input and construction mode.
        if (!this.parentDictionary.features.containsKey(this.featureIdentfier)) {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Feature "
                            + this.featureIdentfier
                            + " not stored in specified parent dictionary"
            );
        }
        if (this.mode.equals(MODE_AMINOACID)
                && !this.parentDictionary.features.get(this.featureIdentfier).isCodingSequence) {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Non-coding feature "
                            + this.featureIdentfier
                            + " incompatible with mode "
                            + this.mode);
        }
        Logging.logStatus(
                "Construct VariantsTable of feature "
                        + featureIdentfier
                        + " and "
                        + sampleIdentifiers.size()
                        + " samples in mode "
                        + this.mode
        );
        switch (mode) {
            case MODE_NUCLEOTIDE -> constructNucleotideTable();
            case MODE_AMINOACID -> constructAminoacidTable();
            default -> throw new MusialException("(VariantsTable) Unknown construction mode " + mode);
        }
    }

    private void constructNucleotideTable() {
        HashMap<Integer, String> perAlleleVariants;
        HashSet<String> processedAlleles = new HashSet<>();
        String sampleAlleleIdentifier;
        int variantStart;
        String variantPosition;
        char[] variantContents;
        char[] referenceContents;
        String variantContent;
        boolean isInsertion;
        boolean isDeletion;
        for (String sampleIdentifier : this.sampleIdentifiers) {
            sampleAlleleIdentifier =
                    this.parentDictionary
                            .samples.get(sampleIdentifier)
                            .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentfier);
            if (sampleAlleleIdentifier.equals(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                continue;
            }
            if (!processedAlleles.contains(sampleAlleleIdentifier)) {
                perAlleleVariants = this.parentDictionary
                        .getNucleotideVariants(this.featureIdentfier, sampleIdentifier, false);
                for (Map.Entry<Integer, String> variantEntry : perAlleleVariants.entrySet()) {
                    variantStart = variantEntry.getKey();
                    variantContents = variantEntry.getValue().toCharArray();
                    referenceContents = this.parentDictionary
                            .nucleotideVariants.get(variantStart).get(variantEntry.getValue())
                            .annotations.get(NucleotideVariantEntry.PROPERTY_NAME_REFERENCE_CONTENT).toCharArray();
                    isInsertion = variantContents.length > referenceContents.length;
                    isDeletion = variantEntry.getValue().contains(String.valueOf(Bio.DELETION_AA1));
                    if ((isInsertion || isDeletion) && this.excludeIndels) {
                        continue;
                    }
                    if (isInsertion) {
                        referenceContents = variantContents;
                        for (int i = 1; i < referenceContents.length; i++) {
                            referenceContents[i] = Bio.DELETION_AA1;
                        }
                    }
                    for (int i = 0; i < variantContents.length; i++) {
                        variantContent = String.valueOf(variantContents[i]);
                        if (isInsertion && i > 0) {
                            variantPosition = variantStart + "+" + i;
                        } else {
                            variantPosition = variantStart + i + "+0";
                        }
                        if (!this.variantsTable.containsKey(variantPosition)) {
                            this.variantsTable.put(variantPosition, new HashMap<>());
                        }
                        if (!this.variantsTable
                                .get(variantPosition)
                                .containsKey(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                            this.variantsTable
                                    .get(variantPosition)
                                    .put(AlleleEntry.PROPERTY_NAME_REFERENCE_ID, String.valueOf(referenceContents[i]));
                        }
                        this.variantsTable
                                .get(variantPosition)
                                .put(sampleAlleleIdentifier, variantContent);
                    }
                }
                processedAlleles.add(sampleAlleleIdentifier);
            }
        }
        // Impute reference content for variant positions per allele.
        for (String processedVariantPosition : this.variantsTable.keySet()) {
            for (String processedAllele : processedAlleles) {
                if (!this.variantsTable
                        .get(processedVariantPosition)
                        .containsKey(processedAllele)) {
                    this.variantsTable
                            .get(processedVariantPosition)
                            .put(
                                    processedAllele,
                                    redundantReferenceContent ?
                                            this.variantsTable.get(processedVariantPosition).get(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)
                                            : "."
                            );
                }
            }
        }
        // Impute reference content for non-variant positions per allele.
        if (this.includeNonVariantPositions) {
            char[] referenceSequenceChars = this.parentDictionary
                    .features.get(this.featureIdentfier)
                    .nucleotideSequence.toCharArray();
            String position;
            String referenceContent;
            for (int i = 0; i < referenceSequenceChars.length; i++) {
                referenceContent = String.valueOf(referenceSequenceChars[i]);
                position = this.parentDictionary.features.get(this.featureIdentfier).start + i + 1 + "+0";
                if (!this.variantsTable.containsKey(position)) {
                    this.variantsTable.put(position, new HashMap<>());
                }
                if (!this.variantsTable
                        .get(position)
                        .containsKey(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                    this.variantsTable
                            .get(position)
                            .put(AlleleEntry.PROPERTY_NAME_REFERENCE_ID, referenceContent);
                }
                for (String processedAllele : processedAlleles) {
                    if (!variantsTable
                            .get(position)
                            .containsKey(processedAllele)) {
                        this.variantsTable
                                .get(position)
                                .put(processedAllele, redundantReferenceContent ? referenceContent : ".");
                    }
                }
            }
        }
    }

    private void constructAminoacidTable() {
        HashMap<String, String> perProteoformVariants;
        HashSet<String> processedProteoforms = new HashSet<>();
        String sampleProteoformIdentifier;
        String variantPosition;
        String variantContent;
        String referenceContent;
        boolean isInsertion;
        boolean isDeletion;
        for (String sampleIdentifier : this.sampleIdentifiers) {
            sampleProteoformIdentifier = this.parentDictionary
                    .samples.get(sampleIdentifier)
                    .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentfier);
            if (sampleProteoformIdentifier.equals(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID)) {
                continue;
            }
            if (!processedProteoforms.contains(sampleProteoformIdentifier)) {
                perProteoformVariants = this.parentDictionary
                        .getProteoformAminoacidVariants(this.featureIdentfier, sampleIdentifier);
                for (Map.Entry<String, String> variantEntry : perProteoformVariants.entrySet()) {
                    variantPosition = variantEntry.getKey();
                    variantContent = variantEntry.getValue();
                    referenceContent = this.parentDictionary
                            .features.get(this.featureIdentfier)
                            .aminoacidVariants.get(variantPosition).get(variantContent)
                            .annotations.get(AminoacidVariantEntry.PROPERTY_NAME_REFERENCE_CONTENT);
                    isInsertion = !variantPosition.split("\\+")[1].equals("0");
                    isDeletion = variantContent.equals(String.valueOf(Bio.DELETION_AA1));
                    if ((isInsertion || isDeletion) && this.excludeIndels) {
                        continue;
                    }
                    if (!this.variantsTable.containsKey(variantPosition)) {
                        this.variantsTable.put(variantPosition, new HashMap<>());
                    }
                    if (!this.variantsTable
                            .get(variantPosition)
                            .containsKey(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                        this.variantsTable
                                .get(variantPosition)
                                .put(AlleleEntry.PROPERTY_NAME_REFERENCE_ID, referenceContent);
                    }
                    this.variantsTable
                            .get(variantPosition)
                            .put(sampleProteoformIdentifier, variantContent);
                }
                processedProteoforms.add(sampleProteoformIdentifier);
            }
        }
        // Impute reference content for variant positions per allele.
        for (String processedVariantPosition : this.variantsTable.keySet()) {
            for (String processedProteoform : processedProteoforms) {
                if (!this.variantsTable
                        .get(processedVariantPosition)
                        .containsKey(processedProteoform)) {
                    this.variantsTable
                            .get(processedVariantPosition)
                            .put(
                                    processedProteoform,
                                    redundantReferenceContent ?
                                            this.variantsTable.get(processedVariantPosition).get(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID)
                                            : "."
                            );
                }
            }
        }
        // Impute reference content for non-variant positions per allele.
        if (this.includeNonVariantPositions) {
            char[] referenceSequenceChars = this.parentDictionary
                    .features.get(this.featureIdentfier)
                    .translatedNucleotideSequence.toCharArray();
            String position;
            for (int i = 0; i < referenceSequenceChars.length; i++) {
                referenceContent = String.valueOf(referenceSequenceChars[i]);
                position = i + 1 + "+0";
                if (!this.variantsTable.containsKey(position)) {
                    this.variantsTable.put(position, new HashMap<>());
                }
                if (!this.variantsTable
                        .get(position)
                        .containsKey(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID)) {
                    this.variantsTable
                            .get(position)
                            .put(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID, referenceContent);
                }
                for (String processedProteoform : processedProteoforms) {
                    if (!variantsTable
                            .get(position)
                            .containsKey(processedProteoform)) {
                        this.variantsTable
                                .get(position)
                                .put(processedProteoform, redundantReferenceContent ? referenceContent : ".");
                    }
                }
            }
        }

    }

    public String toString(boolean grouped) throws MusialException {
        if (!this.mode.equals(MODE_AMINOACID) && !this.mode.equals(MODE_NUCLEOTIDE)) {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Unable to convert table to string with mode "
                            + this.mode
            );
        }
        StringBuilder tableStringBuilder = new StringBuilder();
        ArrayList<String> fields = new ArrayList<>();
        fields.add("Position");
        if (grouped) {
            HashSet<String> groupIdentifiers = new HashSet<>();
            this.sampleIdentifiers.forEach(sId -> groupIdentifiers.add(
                    switch (this.mode) {
                        case MODE_NUCLEOTIDE -> this.parentDictionary
                                .samples.get(sId)
                                .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentfier);
                        case MODE_AMINOACID -> this.parentDictionary
                                .samples.get(sId)
                                .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentfier);
                    }
            ));
            fields.add(
                    switch (this.mode) {
                        case MODE_NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                        case MODE_AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                    }
            );
            fields.addAll(groupIdentifiers);
            tableStringBuilder
                    .append(String.join("\t", fields))
                    .append("\n");
            for (String storedPosition : this.variantsTable.keySet()) {
                fields.clear();
                fields.add(storedPosition);
                fields.add(
                        this.variantsTable.get(storedPosition).get(
                                switch (this.mode) {
                                    case MODE_NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                                    case MODE_AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                                }
                        )
                );
                for (String groupIdentifier : groupIdentifiers) {
                    fields.add(
                            this.variantsTable.get(storedPosition).get(groupIdentifier)
                    );
                }
                tableStringBuilder
                        .append(String.join("\t", fields))
                        .append("\n");
            }
        } else {
            fields.add("Reference");
            fields.addAll(this.sampleIdentifiers);
            tableStringBuilder
                    .append(String.join("\t", fields))
                    .append("\n");
            for (String storedPosition : this.variantsTable.keySet()) {
                fields.clear();
                fields.add(storedPosition);
                fields.add(
                        this.variantsTable.get(storedPosition).get(
                                switch (this.mode) {
                                    case MODE_NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                                    case MODE_AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                                }
                        )
                );
                for (String sampleIdentifier : this.sampleIdentifiers) {
                    fields.add(
                            this.variantsTable.get(storedPosition).get(
                                    switch (this.mode) {
                                        case MODE_NUCLEOTIDE -> this.parentDictionary
                                                .samples.get(sampleIdentifier)
                                                .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentfier);
                                        case MODE_AMINOACID -> this.parentDictionary
                                                .samples.get(sampleIdentifier)
                                                .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentfier);
                                    }
                            )
                    );
                }
                tableStringBuilder
                        .append(String.join("\t", fields))
                        .append("\n");
            }
        }
        return tableStringBuilder.toString();
    }


    @SuppressWarnings("UnusedAssignment")
    public Tuple<String, String> getFastaEntry(String accessIdentifier) throws MusialException {
        String entryIdentifier;
        boolean isSample;
        if (!this.mode.equals(MODE_AMINOACID) && !this.mode.equals(MODE_NUCLEOTIDE)) {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Unable to generate .fasta entry with mode "
                            + this.mode
            );
        }
        if (this.parentDictionary.samples.containsKey(accessIdentifier)) {
            entryIdentifier = switch (this.mode) {
                case MODE_NUCLEOTIDE -> this.parentDictionary.samples.get(accessIdentifier)
                        .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentfier);
                case MODE_AMINOACID -> this.parentDictionary.samples.get(accessIdentifier)
                        .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentfier);
            };
            isSample = true;
        } else if (
                (this.mode.equals(MODE_NUCLEOTIDE) && this.parentDictionary.features.get(this.featureIdentfier).alleles.containsKey(accessIdentifier))
                        || (this.mode.equals(MODE_AMINOACID) && this.parentDictionary.features.get(this.featureIdentfier).proteoforms.containsKey(accessIdentifier))) {
            entryIdentifier = accessIdentifier;
            isSample = false;
        } else {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Unable to obtain .fasta entry for accessor "
                            + accessIdentifier
            );
        }
        String header = ">lcl|"
                + accessIdentifier
                + " [FEATURE=" + this.featureIdentfier + "]";
        String sequence = this.variantsTable.values().stream().map(hm -> hm.get(entryIdentifier)).collect(Collectors.joining());
        return new Tuple<>(header, sequence);
    }
}
