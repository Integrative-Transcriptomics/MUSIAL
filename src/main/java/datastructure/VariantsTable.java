package datastructure;

import cli.ModuleExtractContentModes;
import components.Bio;
import components.Logging;
import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;

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

    public final boolean redundantReferenceContent;

    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(Variants Table)";

    public VariantsTable(VariantsDictionary parentDictionary, ArrayList<String> sampleIdentifiers,
                         String featureIdentifier, ModuleExtractContentModes contentMode, boolean excludeIndels,
                         boolean includeNonVariantPositions, boolean redundantReferenceContent) throws MusialException {
        this.parentDictionary = parentDictionary;
        this.sampleIdentifiers = this.parentDictionary.removeInvalidSampleIds(sampleIdentifiers);
        this.featureIdentifier = featureIdentifier;
        this.contentMode = contentMode;
        this.excludeIndels = excludeIndels;
        this.includeNonVariantPositions = includeNonVariantPositions;
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
        if (this.contentMode.equals(ModuleExtractContentModes.AMINOACID)
                && !this.parentDictionary.features.get(this.featureIdentifier).considerCodingSequence) {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Non-coding feature "
                            + this.featureIdentifier
                            + " incompatible with mode "
                            + this.contentMode);
        }
        Logging.logStatus(
                "Construct variants table ("
                        + this.contentMode
                        + ") of feature "
                        + featureIdentifier
                        + " and "
                        + sampleIdentifiers.size()
                        + " samples"
        );
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

    public String toString(boolean grouped) {
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
                                switch (this.contentMode) {
                                    case NUCLEOTIDE -> AlleleEntry.PROPERTY_NAME_REFERENCE_ID;
                                    case AMINOACID -> ProteoformEntry.PROPERTY_NAME_REFERENCE_ID;
                                }
                        )
                );
                for (String sampleIdentifier : this.sampleIdentifiers) {
                    fields.add(
                            this.variantsTable.get(storedPosition).get(
                                    switch (this.contentMode) {
                                        case NUCLEOTIDE -> this.parentDictionary
                                                .samples.get(sampleIdentifier)
                                                .annotations.get("AL" + FIELD_SEPARATOR_1 + this.featureIdentifier);
                                        case AMINOACID -> this.parentDictionary
                                                .samples.get(sampleIdentifier)
                                                .annotations.get("PF" + FIELD_SEPARATOR_1 + this.featureIdentifier);
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

    public Tuple<String, String> getFastaEntry(String accessIdentifier) throws MusialException {
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
        } else {
            throw new MusialException(
                    EXCEPTION_PREFIX
                            + " Unable to obtain .fasta entry for accessor "
                            + accessIdentifier
            );
        }
        StringBuilder headerStringBuilder = new StringBuilder();
        headerStringBuilder
                .append(">lcl|")
                .append(accessIdentifier)
                .append(" [FEATURE=")
                .append(this.featureIdentifier)
                .append("]");
        HashMap<String, String> annotations = switch (this.contentMode) {
            case NUCLEOTIDE -> this.parentDictionary.features.get(featureIdentifier)
                    .alleles.get(accessIdentifier).annotations;
            case AMINOACID -> this.parentDictionary.features.get(featureIdentifier)
                    .proteoforms.get(accessIdentifier).annotations;
        };
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
            sequenceStringBuilder.append(value.get(entryIdentifier));
        }
        return new Tuple<>(headerStringBuilder.toString(), sequenceStringBuilder.toString());
    }
}
