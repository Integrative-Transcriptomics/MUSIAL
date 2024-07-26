package utility;

import datastructure.*;
import main.Constants;
import org.apache.commons.lang3.tuple.Pair;
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
    @SuppressWarnings({"FieldCanBeLocal", "MismatchedQueryAndUpdateOfCollection"})
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
        Iterator<String> iterator = _featureCoding.getAlleleNameIterator();
        sites = new ArrayList<>();
        ArrayList<String> alleleSequenceContainer; // Stores per position content as one or multiple base symbols incl. gaps for deletions.
        TreeMap<Integer, Pair<String, VariantInformation>> alleleVariants;
        String referenceSequence;
        String proteoformSequence;
        String ref;
        String alt;
        String type;
        int pos;
        int maxInDelLength = 0;
        int relativePosition;
        int alignmentBandWidth;
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
                    referenceSequence = musialStorage.getReferenceSequenceOfFeature(_featureCoding.name);
                    alleleVariants = _featureCoding.getNucleotideVariantsOf(allele.name);
                    alleleSequenceContainer = new ArrayList<>(Arrays.asList(
                            referenceSequence.split("")
                    ));
                    for (Map.Entry<Integer, Pair<String, VariantInformation>> alleleVariantsEntry : alleleVariants.entrySet()) {
                        pos = alleleVariantsEntry.getKey();
                        alt = alleleVariantsEntry.getValue().getLeft();
                        ref = alleleVariantsEntry.getValue().getRight().referenceContent;
                        relativePosition = pos - _featureCoding.start; // 0-based position on feature.
                        type = alleleVariantsEntry.getValue().getRight().getInfo(Constants.TYPE);
                        if (relativePosition < 0) {
                            Logger.logWarning("Skipping variant " + alt + " at (relative) position " + pos + " for allele " + allele.name + " of feature " + featureCoding.name + ".");
                            continue;
                        }
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
        ArrayList<String> sitesList = new ArrayList<>(sites.size());
        String sitesString = "";
        String formName;
        if (sites.size() == 0) {
            formName = Constants.REFERENCE_FORM_NAME;
        } else {
            for (Triple<Integer, String, String> s : sites) {
                final int pos = s.getLeft();
                final String ref = s.getMiddle();
                final String alt = s.getRight();
                sitesList.add(pos + Constants.FIELD_SEPARATOR + alt);
                _featureCoding.addAminoacidVariant(pos, alt, ref);
                allele.getOccurrenceIterator().forEachRemaining(sampleName -> _featureCoding.getAminoacidVariant(pos, alt).addOccurrence(sampleName));
            }
            sitesString = String.join(Constants.ENTRY_SEPARATOR, sitesList);
            formName = storage.getFormName(
                    sitesString.hashCode(),
                    "P" + (_featureCoding.getProteoformCount() + 1)
            );
        }
        Form proteoform = new Form(formName, sitesString);
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
