package datastructure;

import cli.ModuleExtractContentModes;
import components.Bio;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.*;

import static datastructure.VariantsDictionary.FIELD_SEPARATOR_1;

public class VariantsTable {

    private final TreeMap<String, HashMap<String, String>> variantsTable = new TreeMap<>((s1, s2) -> {
        int p1 = Integer.parseInt(s1.split("\\+")[0]);
        int p2 = Integer.parseInt(s2.split("\\+")[0]);
        if (p1 != p2) {
            return Integer.compare(p1, p2);
        } else {
            int i1 = Integer.parseInt(s1.split("\\+")[1]);
            int i2 = Integer.parseInt(s2.split("\\+")[1]);
            return Integer.compare(i1, i2);
        }
    });

    private final VariantsDictionary parentDictionary;

    private final ArrayList<String> sampleIdentifiers;

    private final String featureIdentifier;

    public final ModuleExtractContentModes contentMode;

    public final boolean excludeIndels;

    public final boolean includeNonVariantPositions;

    public HashMap<String, ArrayList<String>> keepOnlyVariantsWith;

    public final boolean redundantReferenceContent;

    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(Variants Table)";

    public VariantsTable(VariantsDictionary parentDictionary, ArrayList<String> sampleIdentifiers,
                         String featureIdentifier, ModuleExtractContentModes contentMode, boolean excludeIndels,
                         boolean includeNonVariantPositions, HashMap<String, ArrayList<String>> keepOnlyVariantsWith, boolean redundantReferenceContent) throws MusialException {
        this.parentDictionary = parentDictionary;
        this.sampleIdentifiers = this.parentDictionary.removeInvalidSampleIds(sampleIdentifiers);
        this.featureIdentifier = featureIdentifier;
        this.contentMode = contentMode;
        this.excludeIndels = excludeIndels;
        this.includeNonVariantPositions = includeNonVariantPositions;
        this.keepOnlyVariantsWith = keepOnlyVariantsWith;
        this.redundantReferenceContent = redundantReferenceContent;
        // Check validity of specified input and construction mode.
        if (!this.parentDictionary.features.containsKey(this.featureIdentifier)) {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Feature "
                            + this.featureIdentifier
                            + " not stored in specified parent dictionary"
            );
        }
        switch (contentMode) {
            case NUCLEOTIDE -> constructNucleotideTable();
            case AMINOACID -> constructAminoacidTable();
        }
    }

    private void constructNucleotideTable() throws MusialException {
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
                            .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentifier);
            if (sampleAlleleIdentifier.equals(AlleleEntry.PROPERTY_NAME_REFERENCE_ID)) {
                continue;
            }
            if (!processedAlleles.contains(sampleAlleleIdentifier)) {
                perAlleleVariants = this.parentDictionary
                        .getNucleotideVariants(this.featureIdentifier, sampleIdentifier, false);
                // Filter variants based on specified annotation values (if any).
                HashSet<Integer> positionsToRemove = new HashSet<>();
                for (Map.Entry<Integer, String> alleleVariantEntry : perAlleleVariants.entrySet()) {
                    HashMap<String, String> alleleVariantAnnotations = parentDictionary.nucleotideVariants.get(alleleVariantEntry.getKey()).get(alleleVariantEntry.getValue()).annotations;
                    // If any variant annotation value was not specified by the user, remove the variant.
                    //noinspection DuplicatedCode
                    for (Map.Entry<String, String> alleleVariantAnnotationEntry : alleleVariantAnnotations.entrySet()) {
                        if (keepOnlyVariantsWith.containsKey(alleleVariantAnnotationEntry.getKey())) {
                            if (!keepOnlyVariantsWith.get(alleleVariantAnnotationEntry.getKey()).contains(alleleVariantAnnotationEntry.getValue())) {
                                positionsToRemove.add(alleleVariantEntry.getKey());
                                break;
                            }
                        }
                    }
                }
                positionsToRemove.forEach(perAlleleVariants::remove);
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
                        referenceContents = new char[variantContents.length];
                        for (int i = 0; i < variantContents.length; i++) {
                            if (i == 0) {
                                referenceContents[i] = variantContents[i];
                            } else {
                                referenceContents[i] = Bio.DELETION_AA1;
                            }
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
        // Impute reference content for non-variant positions per allele.
        if (this.includeNonVariantPositions) {
            char[] referenceSequenceChars = this.parentDictionary
                    .features.get(this.featureIdentifier)
                    .nucleotideSequence.toCharArray();
            String position;
            String referenceContent;
            for (int i = 0; i < referenceSequenceChars.length; i++) {
                referenceContent = String.valueOf(referenceSequenceChars[i]);
                position = this.parentDictionary.features.get(this.featureIdentifier).start + i + "+0";
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
            }
        }
    }

    private void constructAminoacidTable() throws MusialException {
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
                    .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentifier);
            if (sampleProteoformIdentifier.equals(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID)) {
                continue;
            }
            if (!processedProteoforms.contains(sampleProteoformIdentifier)) {
                perProteoformVariants = this.parentDictionary
                        .getAminoacidVariants(this.featureIdentifier, sampleIdentifier);
                // Filter variants based on specified annotation values (if any).
                HashSet<String> positionsToRemove = new HashSet<>();
                for (Map.Entry<String, String> proteoformVariantEntry : perProteoformVariants.entrySet()) {
                    HashMap<String, String> proteoformVariantAnnotations = parentDictionary.features.get(featureIdentifier).aminoacidVariants.get(proteoformVariantEntry.getKey()).get(proteoformVariantEntry.getValue()).annotations;
                    // If any variant annotation value was not specified by the user, remove the variant.
                    //noinspection DuplicatedCode
                    for (Map.Entry<String, String> proteoformVariantAnnotationsEntry : proteoformVariantAnnotations.entrySet()) {
                        if (keepOnlyVariantsWith.containsKey(proteoformVariantAnnotationsEntry.getKey())) {
                            if (!keepOnlyVariantsWith.get(proteoformVariantAnnotationsEntry.getKey()).contains(proteoformVariantAnnotationsEntry.getValue())) {
                                positionsToRemove.add(proteoformVariantEntry.getKey());
                                break;
                            }
                        }
                    }
                }
                positionsToRemove.forEach(perProteoformVariants::remove);
                for (Map.Entry<String, String> variantEntry : perProteoformVariants.entrySet()) {
                    variantPosition = variantEntry.getKey();
                    variantContent = variantEntry.getValue();
                    referenceContent = this.parentDictionary
                            .features.get(this.featureIdentifier)
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
                            .containsKey(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID)) {
                        this.variantsTable
                                .get(variantPosition)
                                .put(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID, referenceContent);
                    }
                    this.variantsTable
                            .get(variantPosition)
                            .put(sampleProteoformIdentifier, variantContent);
                }
                processedProteoforms.add(sampleProteoformIdentifier);
            }
        }
        // Impute reference content for non-variant positions per allele.
        if (this.includeNonVariantPositions) {
            char[] referenceSequenceChars = this.parentDictionary
                    .features.get(this.featureIdentifier)
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
            }
        }
    }

    public void writeTableToFile(String fileName, boolean grouped) throws IOException {
        try (Writer outputWriter = new BufferedWriter(new FileWriter(fileName))) {
            StringBuilder tableStringBuilder = new StringBuilder();
            ArrayList<String> fields = new ArrayList<>();
            fields.add("Position");
            if (grouped) {
                HashSet<String> groupIdentifiers = new HashSet<>();
                this.sampleIdentifiers.forEach(sId -> groupIdentifiers.add(
                        switch (this.contentMode) {
                            case NUCLEOTIDE -> this.parentDictionary
                                    .samples.get(sId)
                                    .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                            case AMINOACID -> this.parentDictionary
                                    .samples.get(sId)
                                    .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                        }
                ));
                groupIdentifiers.remove(
                        switch (this.contentMode) {
                            case NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                            case AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                        }
                );
                fields.add(
                        switch (this.contentMode) {
                            case NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                            case AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                        }
                );
                fields.addAll(groupIdentifiers);
                tableStringBuilder
                        .append(String.join("\t", fields))
                        .append("\n");
                outputWriter.write(tableStringBuilder.toString());
                tableStringBuilder.setLength(0);
                for (String storedPosition : this.variantsTable.keySet()) {
                    fields.clear();
                    fields.add(storedPosition);
                    fields.add(
                            this.variantsTable.get(storedPosition).get(
                                    switch (this.contentMode) {
                                        case NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                                        case AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                                    }
                            )
                    );
                    for (String groupIdentifier : groupIdentifiers) {
                        fields.add(
                                this.variantsTable.get(storedPosition).getOrDefault(groupIdentifier, ".")
                        );
                    }
                    tableStringBuilder
                            .append(String.join("\t", fields))
                            .append("\n");
                    outputWriter.write(tableStringBuilder.toString());
                    tableStringBuilder.setLength(0);
                }
            } else {
                fields.add("Reference");
                fields.addAll(this.sampleIdentifiers);
                tableStringBuilder
                        .append(String.join("\t", fields))
                        .append("\n");
                outputWriter.write(tableStringBuilder.toString());
                tableStringBuilder.setLength(0);
                String sampleGroupIdentifier;
                for (String storedPosition : this.variantsTable.keySet()) {
                    fields.clear();
                    fields.add(storedPosition);
                    fields.add(
                            this.variantsTable.get(storedPosition).get(
                                    switch (this.contentMode) {
                                        case NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                                        case AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                                    }
                            )
                    );
                    for (String sampleIdentifier : this.sampleIdentifiers) {
                        sampleGroupIdentifier = switch (this.contentMode) {
                            case NUCLEOTIDE -> this.parentDictionary
                                    .samples.get(sampleIdentifier)
                                    .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                            case AMINOACID -> this.parentDictionary
                                    .samples.get(sampleIdentifier)
                                    .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                        };
                        if (!Objects.equals(sampleGroupIdentifier, switch (this.contentMode) {
                            case NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                            case AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                        })) {
                            fields.add(
                                    this.variantsTable.get(storedPosition).getOrDefault(
                                            switch (this.contentMode) {
                                                case NUCLEOTIDE -> this.parentDictionary
                                                        .samples.get(sampleIdentifier)
                                                        .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                                                case AMINOACID -> this.parentDictionary
                                                        .samples.get(sampleIdentifier)
                                                        .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                                            },
                                            "."
                                    )
                            );
                        } else {
                            fields.add(".");
                        }
                    }
                    tableStringBuilder
                            .append(String.join("\t", fields))
                            .append("\n");
                    outputWriter.write(tableStringBuilder.toString());
                    tableStringBuilder.setLength(0);
                }
            }
        }
    }

    public Tuple<String, String> getFastaEntry(String accessIdentifier, boolean removeGaps) throws MusialException {
        boolean accessIsReference = (accessIdentifier.equals(AlleleEntry.PROPERTY_NAME_REFERENCE_ID) || accessIdentifier.equals(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID));
        boolean referenceIsStored = (this.contentMode.equals(ModuleExtractContentModes.NUCLEOTIDE) && this.parentDictionary.features.get(this.featureIdentifier).alleles.containsKey(AlleleEntry.PROPERTY_NAME_REFERENCE_ID))
                || (this.contentMode.equals(ModuleExtractContentModes.AMINOACID) && this.parentDictionary.features.get(this.featureIdentifier).proteoforms.containsKey(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID));
        String entryIdentifier;
        if (this.parentDictionary.samples.containsKey(accessIdentifier)) {
            entryIdentifier = switch (this.contentMode) {
                case NUCLEOTIDE -> this.parentDictionary.samples.get(accessIdentifier)
                        .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                case AMINOACID -> this.parentDictionary.samples.get(accessIdentifier)
                        .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentifier);
            };
        } else if (
                (this.contentMode.equals(ModuleExtractContentModes.NUCLEOTIDE) && this.parentDictionary.features.get(this.featureIdentifier).alleles.containsKey(accessIdentifier))
                        || (this.contentMode.equals(ModuleExtractContentModes.AMINOACID) && this.parentDictionary.features.get(this.featureIdentifier).proteoforms.containsKey(accessIdentifier))) {
            entryIdentifier = accessIdentifier;
        } else if (!accessIsReference) {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Unable to obtain .fasta entry for accessor "
                            + accessIdentifier
            );
        } else {
            entryIdentifier = accessIdentifier;
        }
        StringBuilder headerStringBuilder = new StringBuilder();
        headerStringBuilder
                .append(">lcl|")
                .append(accessIdentifier)
                .append(" [FEATURE=")
                .append(this.featureIdentifier)
                .append("]");
        HashMap<String, String> annotations;
        if (accessIsReference && !referenceIsStored) {
            annotations = new HashMap<>() {{
                put(AlleleEntry.PROPERTY_NAME_FREQUENCY, "0.00");
            }};
        } else if (accessIsReference) {
            annotations = new HashMap<>() {{
                put(
                        AlleleEntry.PROPERTY_NAME_FREQUENCY,
                        switch (contentMode) {
                            case NUCLEOTIDE -> parentDictionary.features.get(featureIdentifier)
                                    .alleles.get(AlleleEntry.PROPERTY_NAME_REFERENCE_ID).annotations.get(AlleleEntry.PROPERTY_NAME_FREQUENCY);
                            case AMINOACID -> parentDictionary.features.get(featureIdentifier)
                                    .proteoforms.get(ProteoformEntry.PROPERTY_NAME_REFERENCE_ID).annotations.get(ProteoformEntry.PROPERTY_NAME_FREQUENCY);
                        }
                );
            }};
        } else {
            annotations = switch (this.contentMode) {
                case NUCLEOTIDE -> this.parentDictionary.features.get(featureIdentifier)
                        .alleles.get(entryIdentifier).annotations;
                case AMINOACID -> this.parentDictionary.features.get(featureIdentifier)
                        .proteoforms.get(entryIdentifier).annotations;
            };
        }
        for (Map.Entry<String, String> annotation : annotations.entrySet()) {
            //noinspection ConstantConditions
            if (annotation.getKey().equals(AlleleEntry.PROPERTY_NAME_VARIANTS)
                    || annotation.getKey().equals(ProteoformEntry.PROPERTY_NAME_VARIANTS)) {
                continue;
            }
            headerStringBuilder
                    .append(" [")
                    .append(annotation.getKey())
                    .append("=")
                    .append(annotation.getValue())
                    .append("]");
        }
        StringBuilder sequenceStringBuilder = new StringBuilder();
        for (HashMap<String, String> value : variantsTable.values()) {
            if (value.containsKey(entryIdentifier)) {
                sequenceStringBuilder.append(value.get(entryIdentifier));
            } else {
                sequenceStringBuilder.append(value.get(switch (this.contentMode) {
                    case NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                    case AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                }));
            }
        }
        String sequenceString = sequenceStringBuilder.toString();
        return new Tuple<>(headerStringBuilder.toString(), removeGaps ? sequenceString.replace("-", "") : sequenceString);
    }
}