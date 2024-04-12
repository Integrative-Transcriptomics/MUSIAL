package utility;

import datastructure.*;
import main.Constants;
import org.apache.commons.lang3.tuple.Triple;
import org.javatuples.Triplet;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableSet;

/**
 * Implements static methods for allele analysis to infer proteoforms.
 * <p>
 * Analyzes alleles of features, aligns them with reference sequences, and derives variant information.
 *
 * @author Simon Hackl
 * @version 2.3
 * @since 2.2
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
        StringBuilder alleleSequenceBuilder;
        String proteoformSequence;
        String ref;
        String alt;
        boolean isDeletion;
        boolean isInsertion;
        boolean isInDel;
        long noDeletedPositions;
        int lengthDifference;
        int relativeVariantPosition;
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
                    alleleSequenceBuilder = new StringBuilder(musialStorage.getReferenceSequenceOfFeature(featureCoding.name));
                    variantAlleleSites = featureCoding.getNucleotideVariantPositions();
                    for (int pos : variantAlleleSites) {
                        for (Map.Entry<String, VariantInformation> entry : featureCoding.getNucleotideVariantsAt(pos).entrySet()) {
                            if (allele.getOccurrenceSet().stream().noneMatch(entry.getValue()::hasOccurrence))
                                continue;
                            relativeVariantPosition = pos - featureCoding.start; // 0-based position on feature.
                            alt = entry.getKey();
                            ref = entry.getValue().referenceContent;
                            isDeletion = alt.contains(Constants.DELETION_OR_GAP_STRING);
                            isInsertion = ref.contains(Constants.DELETION_OR_GAP_STRING);
                            isInDel = (isDeletion && isInsertion);
                            lengthDifference = Math.abs(ref.replaceAll("-", "").length() - alt.replaceAll("-", "").length());
                            if (isInDel) {
                                // Insert Insertion-Deletion.
                                // Compare with deletion and insertion case; Replaces deleted reference content with put. inserted alt. content.
                                noDeletedPositions = alt.chars().filter(c -> c == '-').count();
                                alleleSequenceBuilder.replace(relativeVariantPosition, (int) (relativeVariantPosition + noDeletedPositions + 1), alt.replace("-", ""));
                                alignmentBandWidth += lengthDifference;
                            } else if (isDeletion) {
                                // Replace deletion.
                                // Replace substitutes substring at [ start, end ).
                                alleleSequenceBuilder.replace(relativeVariantPosition, relativeVariantPosition + alt.length(), alt.replace("-", ""));
                                alignmentBandWidth += lengthDifference;
                            } else if (isInsertion) {
                                // Insert insertion.
                                // Insert adds substring starting at offset, i.e., we need to shift by one to account for excluded conserved position.
                                alleleSequenceBuilder.insert(relativeVariantPosition + 1, alt.substring(1));
                                alignmentBandWidth += lengthDifference;
                            } else {
                                // Case; Replace substitution; shift end by one for fixed subst. length of one.
                                alleleSequenceBuilder.replace(relativeVariantPosition, relativeVariantPosition + alt.length(), alt);
                            }
                        }
                    }
                    alignmentBandWidth = (int) Math.ceil(alignmentBandWidth / 3.0);
                    // Translate allele nucleotide sequence and align with reference aminoacid sequence.
                    proteoformSequence = SequenceOperations.translateNucSequence(alleleSequenceBuilder.toString(), true, true, featureCoding.isSense);
                    Triplet<Integer, String, String> sequenceAlignment = SequenceOperations.globalAminoAcidSequenceAlignment(
                            featureCoding.getCodingSequence(),
                            proteoformSequence,
                            proteoformSequence.length() * 8, 9,
                            SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID,
                            SequenceOperations.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE,
                            alignmentBandWidth
                    );
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
