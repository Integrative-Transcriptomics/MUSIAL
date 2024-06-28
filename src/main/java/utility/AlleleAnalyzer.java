package utility;

import datastructure.*;
import main.Constants;
import org.apache.commons.lang3.tuple.Triple;
import org.javatuples.Triplet;

import java.util.*;

/**
 * Implements static methods for allele analysis to infer proteoforms.
 * <p>
 * Analyzes alleles of features, aligns them with reference sequences, and derives variant information.
 *
 * @author Simon Hackl
 */
public class AlleleAnalyzer {

    /**
     * The MusialStorage instance containing genetic variant data.
     */
    private static MusialStorage storage;

    /**
     * The current FeatureCoding instance to be analyzed.
     */
    private static FeatureCoding _featureCoding;

    /**
     * The current allele being analyzed.
     */
    private static Form allele;

    /**
     * Container of all called sites to infer alleles.
     */
    private static ArrayList<Triple<Integer, String, String>> sites;

    /**
     * List of positions where novel terminations are found.
     */
    private static ArrayList<String> novelTerminations;

    /**
     * Runs the allele analysis for a given FeatureCoding and MusialStorage.
     *
     * @param featureCoding The FeatureCoding instance representing the feature to be analyzed.
     * @param musialStorage The MusialStorage instance containing genetic variant data.
     */
    public static void run(FeatureCoding featureCoding, MusialStorage musialStorage) {
        storage = musialStorage;
        _featureCoding = featureCoding;
        Iterator<String> iterator = featureCoding.getAlleleNameIterator();
        sites = new ArrayList<>();
        NavigableSet<Integer> variantAlleleSites;
        ArrayList<String> alleleSequenceContainer; // Stores per position content as one or multiple base symbols incl. gaps for deletions.
        String referenceSequence;
        String proteoformSequence;
        String ref;
        String alt;
        String type;
        int maxInDelLength = 0;
        int relativePosition;
        int alignmentBandWidth = 1;
        char[] proteoformSequenceAlnChars;
        char proteoformSequenceAlnChar;
        while (iterator.hasNext()) {
            allele = _featureCoding.getAllele(iterator.next());
            try {
                if (allele.name.equals(Constants.REFERENCE_FORM_NAME)) {
                    // Proteoform is equal to the reference.
                    novelTerminations = new ArrayList<>();
                    sites = new ArrayList<>();
                } else {
                    // Derive allele nucleotide sequence and maximal indel length.
                    referenceSequence = musialStorage.getReferenceSequenceOfFeature(featureCoding.name);
                    alleleSequenceContainer = new ArrayList<>(Arrays.asList(
                            referenceSequence.split("")
                    ));
                    variantAlleleSites = featureCoding.getNucleotideVariantPositions(); // TODO: Could be improved by method to return VariantInformation entries by allele name.
                    for (int pos : variantAlleleSites) {
                        for (Map.Entry<String, VariantInformation> entry : featureCoding.getNucleotideVariantsAt(pos).entrySet()) {
                            if (allele.getOccurrenceSet().stream().noneMatch(entry.getValue()::hasOccurrence))
                                continue;
                            relativePosition = pos - featureCoding.start; // 0-based position on feature.
                            alt = entry.getKey();
                            ref = entry.getValue().referenceContent;
                            type = entry.getValue().getInfo(Constants.TYPE);
                            if (type.contains(Constants.TYPE_INSERTION) || type.contains(Constants.TYPE_SUBSTITUTION))
                                alleleSequenceContainer.set(relativePosition, alt);
                            else
                                for (int i = 0; i < alt.length(); i++) {
                                    if ((relativePosition + i) < alleleSequenceContainer.size()) // Avoid deletions that exceed the feature length!
                                        alleleSequenceContainer.set(relativePosition + i, String.valueOf(alt.charAt(i)));
                                }
                            maxInDelLength = Math.max(
                                    maxInDelLength,
                                    Math.abs(ref.replaceAll(Constants.DELETION_OR_GAP_STRING, "").length() - alt.replaceAll(Constants.DELETION_OR_GAP_STRING, "").length())
                            );
                        }
                    }

                    alignmentBandWidth = Math.max(12, (int) Math.ceil(maxInDelLength / 3.0));
                    // Translate allele nucleotide sequence and align with reference aminoacid sequence.
                    proteoformSequence = SequenceOperations.translateNucSequence(
                            String.join("", alleleSequenceContainer).replace(Constants.DELETION_OR_GAP_STRING, ""),
                            true,
                            true,
                            featureCoding.isSense
                    );
                    Triplet<Integer, String, String> sequenceAlignment;
                    try {
                        sequenceAlignment = SequenceOperations.globalAminoAcidSequenceAlignment(
                                featureCoding.getCodingSequence(),
                                proteoformSequence,
                                12, 6,
                                SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID,
                                SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE,
                                alignmentBandWidth
                        );
                    } catch (IndexOutOfBoundsException e) {
                        Logger.logWarning("Banded alignment failed for feature " + featureCoding.name + " and allele " + allele.name + "; Re-running analysis as full alignment.");
                        sequenceAlignment = SequenceOperations.globalAminoAcidSequenceAlignment(
                                featureCoding.getCodingSequence(),
                                proteoformSequence,
                                12, 6,
                                SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID,
                                SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE,
                                0
                        );
                    }
                    // Refine aligned allele aminoacid sequence notation to use for storage.
                    proteoformSequenceAlnChars = sequenceAlignment.getValue2().toCharArray();
                    novelTerminations = new ArrayList<>();
                    // Check proteoform sequence for novel terminations.
                    for (int i = 0; i < proteoformSequenceAlnChars.length - 1; i++) {
                        proteoformSequenceAlnChar = proteoformSequenceAlnChars[i];
                        if (proteoformSequenceAlnChar == Constants.TERMINATION_AMINOACID_CHAR) {
                            novelTerminations.add(String.valueOf(i + 1));
                        }
                    }
                    sites = SequenceOperations.getVariantsOfAlignedSequences(sequenceAlignment.getValue1(), sequenceAlignment.getValue2());
                }
                allocateVariantsInformation();
            } catch (Exception e) {
                Logger.logError("An error occurred during the analysis of allele " + allele.name + " of feature " + featureCoding.name + "; " + e.getMessage());
                e.printStackTrace();
            }
        }
    }

    /**
     * Allocates variant information for the analyzed allele.
     */
    @SuppressWarnings("DuplicatedCode")
    private static void allocateVariantsInformation() {
        StringBuilder sitesString = new StringBuilder(sites.size() * 4);
        String formName;
        if (sites.size() == 0) {
            formName = Constants.REFERENCE_FORM_NAME;
        } else {
            int count_substitution = 0;
            int count_insertion = 0;
            int count_deletion = 0;
            int count_ambiguous = 0;
            for (Triple<Integer, String, String> s : sites) {
                final int pos = s.getLeft();
                final String ref = s.getMiddle();
                final String alt = s.getRight();
                if (alt.contains(Constants.ANY_AMINOACID_STRING)) {
                    count_ambiguous += alt.codePoints().filter(c -> c == Constants.ANY_AMINOACID_CHAR).count();
                } else if (ref.length() == 1 && alt.length() == 1) {
                    count_substitution += 1;
                } else {
                    count_insertion += ref.codePoints().filter(c -> c == '-').count();
                    count_deletion += alt.codePoints().filter(c -> c == '-').count();
                }
                sitesString.append(pos).append(Constants.FIELD_SEPARATOR_1).append(alt).append(Constants.FIELD_SEPARATOR_2);
                _featureCoding.addAminoacidVariant(pos, alt, ref);
                allele.getOccurrenceIterator().forEachRemaining(sampleName -> _featureCoding.getAminoacidVariant(pos, alt).addOccurrence(sampleName));
            }
            sitesString.setLength(sitesString.length() - 1); // Delete last ';' separator symbol.
            formName = storage.getFormName(
                    sitesString.toString().hashCode(),
                    "P" + (_featureCoding.getProteoformCount() + 1) + ".s" + count_substitution + ".i" + count_insertion + ".d" + count_deletion + ".a" + count_ambiguous + ".t" + novelTerminations.size()
            );
        }
        Form proteoform = new Form(formName, sitesString.toString());
        _featureCoding.addProteoform(proteoform);
        Iterator<String> iterator = allele.getOccurrenceIterator();
        do {
            Sample sample = storage.getSample(iterator.next());
            _featureCoding.getProteoform(proteoform.name).addOccurrence(sample.name);
            sample.setProteoform(_featureCoding.name, proteoform.name);
        } while (iterator.hasNext());
        allele.addInfo("proteoform", proteoform.name);
    }

}
