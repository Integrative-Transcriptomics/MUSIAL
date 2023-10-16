package runnables;

import datastructure.FeatureCoding;
import datastructure.Form;
import datastructure.MusialStorage;
import main.MusialConstants;
import org.javatuples.Triplet;
import utility.Bio;
import utility.Compression;
import utility.Logger;

import java.util.ArrayList;
import java.util.TreeMap;

/**
 * Implementation of the {@link Runnable} interface to execute threaded allele analysis to infer proteoforms.
 *
 * @author Simon Hackl
 * @version 2.2
 * @since 2.2
 */
public class AlleleAnalyzer implements Runnable {

    /**
     * The {@link FeatureCoding} to be analyzed.
     */
    FeatureCoding featureCoding;

    /**
     * The {@link Form} to be analyzed.
     */
    Form allele;

    /**
     * The {@link MusialStorage} to transfer proteoform information to.
     */
    MusialStorage musialStorage;

    /**
     * Constructor of the {@link AlleleAnalyzer} class.
     *
     * @param featureCoding The {@link FeatureCoding} of which an allele should be analyzed.
     * @param allele        The {@link Form} to be analyzed.
     * @param musialStorage The {@link MusialStorage} to transfer information to.
     */
    public AlleleAnalyzer(FeatureCoding featureCoding, Form allele, MusialStorage musialStorage) {
        this.featureCoding = featureCoding;
        this.allele = allele;
        this.musialStorage = musialStorage;
    }

    @Override
    public void run() {
        try {
            int relativeVariantPosition;
            int maximalIndelLength = 0;
            String variantContent;
            String proteoformName;
            if (allele.name.equals(MusialConstants.REFERENCE_ID)) {
                // Case; The proteoform is equal to the reference.
                proteoformName = MusialConstants.REFERENCE_ID;
                Form proteoform = new Form(proteoformName);
                proteoform.addAnnotation(MusialConstants.VARIANTS, "");
                featureCoding.addProteoform(proteoform);
            } else {
                // Derive allele nucleotide sequence and maximal indel length.
                StringBuilder sequenceBuilder = new StringBuilder(musialStorage.getReferenceSequenceOfFeature(featureCoding.name));
                String referenceContent;
                TreeMap<Integer, String> variants = featureCoding.getNucleotideVariants(allele.name);
                for (Integer variantPosition : variants.navigableKeySet()) {
                    relativeVariantPosition = variantPosition - featureCoding.start; // 0-based position on feature.
                    variantContent = variants.get(variantPosition);
                    referenceContent = featureCoding.getNucleotideVariantAnnotation(variantPosition, variantContent).getProperty(MusialConstants.REFERENCE_CONTENT);
                    boolean isDeletion = variantContent.contains(String.valueOf(Bio.DELETION_AA1));
                    boolean isInsertion = referenceContent.contains(String.valueOf(Bio.DELETION_AA1));
                    boolean isInDel = (isDeletion && isInsertion);
                    if (isInDel) {
                        // Case; Insert InDel.
                        // Compare with deletion and insertion case; Replaces deleted reference content with put. inserted alt. content.
                        long deletions = variantContent.chars().filter(c -> c == '-').count();
                        sequenceBuilder.replace(relativeVariantPosition, (int) (relativeVariantPosition + deletions + 1), variantContent.replace("-", ""));
                    } else if (isDeletion) {
                        // Case; Replace deletion.
                        // Replace substitutes substring at [ start, end ).
                        sequenceBuilder.replace(relativeVariantPosition, relativeVariantPosition + variantContent.length(), variantContent.replace("-", ""));
                    } else if (isInsertion) {
                        // Case; Insert insertion.
                        // Insert adds substring starting at offset, i.e., we need to shift by one to account for excluded conserved position.
                        sequenceBuilder.insert(relativeVariantPosition + 1, variantContent.substring(1));
                    } else {
                        // Case; Replace substitution; shift end by one for fixed subst. length of one.
                        sequenceBuilder.replace(relativeVariantPosition, relativeVariantPosition + 1, variantContent);
                    }
                    maximalIndelLength = Math.max(maximalIndelLength, variantContent.length());
                }
                maximalIndelLength = (int) Math.ceil(maximalIndelLength / 3.0);
                // Translate allele nucleotide sequence and align with reference aminoacid sequence.
                String alleleAminoacidSequence = Bio.translateNucSequence(sequenceBuilder.toString(), true, true, featureCoding.isSense);
                Triplet<Integer, String, String> sequenceAlignment = Bio.globalAminoAcidSequenceAlignment(
                        featureCoding.getCodingSequence(),
                        alleleAminoacidSequence,
                        alleleAminoacidSequence.length() * 8, 11,
                        Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FORBID,
                        Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.PENALIZE,
                        maximalIndelLength
                );
                // Refine aligned allele aminoacid sequence notation to use for storage.
                char[] alignedAlleleAminoacidSequence = sequenceAlignment.getValue2().toCharArray();
                char alignedAlleleAminoacidCharacter;
                boolean hasNovelStop = false;
                // Check proteoform sequence for novel terminations.
                ArrayList<String> novelStops = new ArrayList<>();
                for (int i = 0; i < alignedAlleleAminoacidSequence.length - 1; i++) {
                    alignedAlleleAminoacidCharacter = alignedAlleleAminoacidSequence[i];
                    if (alignedAlleleAminoacidCharacter == Bio.TERMINATION_AA1) {
                        hasNovelStop = true;
                        novelStops.add(String.valueOf(i + 1));
                    }
                }
                ArrayList<String> proteinVariants = Bio.getVariantsOfAlignedSequences(sequenceAlignment.getValue1(), sequenceAlignment.getValue2());
                if (proteinVariants.size() == 0) {
                    // Case; The proteoform is equal to the reference. I.e., allele has nucleotide variants, but they do not change protein sequence.
                    proteoformName = MusialConstants.REFERENCE_ID;
                    Form proteoform = new Form(proteoformName);
                    proteoform.addAnnotation(MusialConstants.VARIANTS, "");
                    featureCoding.addProteoform(proteoform);
                } else {
                    String position;
                    StringBuilder variantFingerprint = new StringBuilder();
                    String[] proteinVariantFields;
                    for (String proteinVariant : proteinVariants) {
                        proteinVariantFields = proteinVariant.split(MusialConstants.FIELD_SEPARATOR_1);
                        position = proteinVariantFields[0];
                        variantContent = proteinVariantFields[1];
                        referenceContent = proteinVariantFields[2];
                        variantFingerprint.append(position).append(MusialConstants.FIELD_SEPARATOR_1).append(variantContent).append(MusialConstants.FIELD_SEPARATOR_2);
                        // Add novel amino acid variant with annotation.
                        featureCoding.addAminoacidVariant(Integer.parseInt(position), variantContent);
                        featureCoding.getAminoacidVariantAnnotation(Integer.parseInt(position), variantContent).addProperty(MusialConstants.REFERENCE_CONTENT, referenceContent);
                        for (String sampleName : featureCoding.getAllele(allele.name).getOccurrence()) {
                            featureCoding.getAminoacidVariantAnnotation(Integer.parseInt(position), variantContent).addProperty(MusialConstants.VARIANT_OCCURRENCE_SAMPLE_PREFIX + sampleName, ""); // FIXME: Any useful content for this?
                        }
                    }
                    if (variantFingerprint.length() != 0) {
                        variantFingerprint.setLength(variantFingerprint.length() - 1);
                    }
                    proteoformName = String.valueOf(variantFingerprint.toString().hashCode());
                    Form proteoform = new Form(proteoformName);
                    proteoform.addAnnotation(MusialConstants.VARIANTS, Compression.brotliEncodeString(variantFingerprint.toString()));
                    if (hasNovelStop)
                        proteoform.addAnnotation(MusialConstants.PROTEOFORM_NOVEL_STOPS, String.join(MusialConstants.FIELD_SEPARATOR_1, novelStops));
                    featureCoding.addProteoform(proteoform);
                }
            }
            // Add allele information.
            featureCoding.getProteoform(proteoformName).addOccurrence(allele.name);
            allele.addAnnotation(MusialConstants.ALLELE_ANNOTATION_PROTEOFORM, proteoformName);
            // Add sample information.
            for (String sampleName : featureCoding.getAllele(allele.name).getOccurrence()) {
                musialStorage.getSample(sampleName).addAnnotation(MusialConstants.SAMPLE_ANNOTATION_PROTEOFORM_PREFIX + featureCoding.name, proteoformName);
            }
        } catch (Exception e) {
            Logger.logError("An error occurred during the analysis of allele " + allele.name + " of feature " + featureCoding.name + "; " + e.getMessage());
            e.printStackTrace();
        }
    }

}
